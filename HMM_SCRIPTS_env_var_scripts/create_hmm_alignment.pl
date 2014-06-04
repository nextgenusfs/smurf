#!/usr/local/bin/perl
##################################################################
# $Header: /usr/local/devel/microbial/src/hmm/hmm_search/src/RCS/create_hmm_alignment.pl,v 1.6 1999/01/14 23:19:56 delwood Exp $
##################################################################
#  Revision History:
#    $Log: create_hmm_alignment.pl,v $
#    Revision 1.6  1999/01/14 23:19:56  delwood
#    modified to allow scores down to -2000.
#
#    Revision 1.5  1998/11/20 16:02:51  delwood
#    revised perl5 to perl, and convert '-' in names to '_' for belvu
#
#    Revision 1.4  1998/08/03 22:16:49  delwood
#    Force names to be shorter for GP names that are too long.
#    THis will allow use of hmmfs hits where coords are important
#
#    Revision 1.3  1998/07/08 01:07:47  delwood
#    added length and score cutoffs to processing of scores.
#    Also did a work around the truncated name problem for GP names.
#
#    Revision 1.2  1998/06/12 18:33:58  delwood
#    Added -M option to save alignment as MSF file.
#    Uses Belvu for conversion.
#
#    Revision 1.1  1998/06/12 14:43:32  delwood
#    Initial revision
#
###############################################################

require "getopts.pl";
$usage = <<_EOT_;
Usage:  create_hmm_alignment.pl -d fasta -c coord_file [ -H hmm_name -s output  ...]
           ____or____
     cat coord_file |  create_hmm_alignment.pl -d fasta [ -H hmm_name -s output ...]

    -d database.fa     mini_db of seqences (required)
    -c coord_file      file of 'name end5 end3' (default: accepts coords on STDIN)
    -s output          output file of trimmed sequences (fasta format, default is tmp file)
    -l length          length of match required
    -b cutoff          bit cutoff for match
    -q                 quiet operation
 ____HMM alignment options____
    -H xxx.HMM	name of hmm to be used for alignment (no alignment if not specifed)
    -O "options" 	hmm alignment options (default is "-q")
    -S align_name       name to save hmm alignment file (default is "database_trim.slx")
    -E                  flag to save score file as "align_name.score"
    -A align_prog       hmm alignment program (default hmmalign)
    -M                  flag to make alignment msf format (from selex)

    It accepts -D for debug mode and -h for help. 
_EOT_


&Getopts('hDd:c:s:l:b:qH:O:S:EA:M');
if($opt_h == 1){
    die "$usage\n";
}
$DEBUG = $opt_D;
if(length($opt_d) > 0){
    $DB = $opt_d;
}else{
    die "Must specify mini database\n$usage\n";
}
if(length($opt_c) > 0){
    $COORD_FILE = $opt_c;
}else{
    $COORD_FILE = "";
}
if(length($opt_s) > 0){
    $SAVE = $opt_s;
}else{
    $SAVE = "/usr/local/projects/family/tmp/$$.db";
}
if(length($opt_l) > 0){
    $LEN_C = $opt_l;
}else{
    $LEN_C = 0;
}
if(length($opt_b) > 0){
    $BIT_C = $opt_b;
}else{
    $BIT_C = -2000;
}
$QUIET = $opt_q;

if(length($opt_H) > 0){
    $HMM = $opt_H;
}else{
    print STDERR "no HMM specified, so not aligning\n";
    $HMM = "";
}
if(length($opt_O) > 0){
    $ALIGN_OPTS = $opt_O;
}else{
    $ALIGN_OPTS = "-q ";
}
if(length($opt_S) > 0){
    $H_SAVE = $opt_S;
}else{
    $H_SAVE = "";
}
if(length($opt_A) > 0){
    $ALIGN_PROG = $opt_A;
}else{
    $ALIGN_PROG = "/usr/local/bin/hmmalign";
}
$MSF_FLAG = $opt_M;

### end options
print STDERR "DEBUG on\n" if $DEBUG;
#first read in mini db
%data = &read_db($DB);

#now read through coords and pull from data
open(D_SAVE, ">$SAVE") || die " can't open save file $SAVE\n";
print STDERR "Saving trimmed sequences in $SAVE\n" if $DEBUG;
if($COORD_FILE ne ""){
    open(F_IN, "$COORD_FILE") || die "can't open cord file $COORD_FILE\n";
    $INPUT = F_IN;
}else{
    $INPUT = STDIN;
}
print STDERR "accepting input from $INPUT\n" if $DEBUG;
while(<$INPUT>){
    chomp;
    @coords = split(/\s+/, $_);
    print STDERR "found name: ..$coords[0]..\n" if $DEBUG;
    if($data{$coords[0]} ne ""){
	$len = $coords[2] - $coords[1] + 1;
	$score = ($coords[3] ne "" ? $coords[3] : 100);
	if($len > $LEN_C && $score >= $BIT_C){
	    $seq = substr $data{$coords[0]}, $coords[1] - 1, $len;
	    ## note that - will cause belvu to cough....
	    $coords[0] =~ s/-/_/g;
	    print D_SAVE ">$coords[0]/$coords[1]-$coords[2]\n";
	    &write_seq(D_SAVE, $seq);
	}else{
	    print STDERR "$coords[0] skipped, score ..$score.. length ..$len..\n" if !$QUIET;
	}
    }else{       # this is for names that were truncated. Output a shorter name here!!
	$real_name = &get_name($coords[0], (keys(%data)));
	print STDERR "Found name ..$real_name.. for $coords[0]\n" if $DEBUG;
	$len = $coords[2] - $coords[1] + 1;
	$score = ($coords[3] ne "" ? $coords[3] : 100);
	if($real_name eq ""){
	    print STDERR "no sequence found for ...$coords[0]... in $DB\n" if !$QUIET;
	}elsif($len > $LEN_C && $score >= $BIT_C){
	    $seq = substr $data{$real_name}, $coords[1] - 1, $len;
	    $short_name = &clip_name($real_name);
	    $short_name =~ s/-/_/g;
	    print D_SAVE ">$short_name/$coords[1]-$coords[2]\n";
	    &write_seq(D_SAVE, $seq);
	}else{
	    print STDERR "$coords[0] skipped, score ..$score.. length ..$len..\n" if !$QUIET;
	}
    }
}
close(F_IN);
close(D_SAVE);

#now align seq to HMM if defined. 
if($HMM ne ""){
    if($HMM =~ /\/?(\w+)\.HMM/  || $HMM =~ /\/?(\w+)\.hmm/){
	$hmm_name = $1;
    }
    if($H_SAVE eq ""){
	$DB =~ s/\.\w+$//;
	if($MSF_FLAG == 1){
	    $H_SAVE = $DB . "_trim.msf";
	    $MSF = "| /usr/local/bin/belvu -o MSF - ";
	}else{
	    $H_SAVE = $DB . "_trim.slx";
	    $MSF = "";
	}

    }elsif($MSF_FLAG == 1){
	$MSF = "| /usr/local/bin/belvu -o MSF - ";
    }
    $SCORE = "-s $H_SAVE.score " if ($opt_E == 1);
    $cmd = "$ALIGN_PROG $ALIGN_OPTS $SCORE $HMM $SAVE $MSF> $H_SAVE";
    system($cmd);
    unlink $SAVE if ($opt_s eq "");
}

exit(0);

#####################################################################
sub read_db{
    local($DB) = @_;
    local(%l, $buf, $id);

    open(DATA, "$DB") || die "can't open $DB for reading.\n";
    while(<DATA>){
	chomp;
	if($_ =~ s/^>//){
	    if($FASTA == 1){
##		print STDERR "name:$id\tseq:$buf\n" if $DEBUG;
		$l{$id} = $buf;
		$id = $buf = "";
	    }else{
		$FASTA = 1;
	    }
	    ($id, $dummy) = split(/\s+/, $_, 2);
	    $id =~ s/\/.+//;
	}elsif($FASTA == 1){
	    $_ =~ s/[\._-]//g;
	    $buf .= $_;
	}
    }
    if($FASTA == 1){
	print STDERR "name:$id\tseq:$buf\n" if $DEBUG;
	$l{$id} = $buf;
    }
    close(DATA);

    return(%l);
}
sub write_seq {
    local($out, $seq)= @_;
    local($i, $len);
    local($FASTA_LEN) = 60;
 
    $len = length ($seq);
    for ($i = 0; $i <= $len; $i += $FASTA_LEN) {
        print $out (substr($seq, $i, $FASTA_LEN)."\n");
    }
}
sub get_name{
    local($name, @list) = @_;
    local(@a);

    @a = grep /\Q$name\E/, @list;
    print STDERR "get_name found ..@a..\n" if $DEBUG;
    return($a[0]);
}

sub clip_name{
    local($name) = @_;
    local(@a, $short_name);

    @a = split(/\|/, $name);
    $short_name = "$a[0]|$a[1]";

    return($short_name);
}
