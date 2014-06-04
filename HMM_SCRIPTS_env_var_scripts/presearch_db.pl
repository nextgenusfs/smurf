#!/usr/local/bin/perl
##################################################################
# $Header: /usr/local/devel/microbial/src/hmm/hmm_search/src/RCS/presearch_db.pl,v 1.7 1999/11/04 16:58:56 delwood Exp $
##################################################################
#  Revision History:
#    $Log: presearch_db.pl,v $
#    Revision 1.7  1999/11/04 16:58:56  delwood
#    revised to use new yank and search against /usr/local/db/panda/nraa/nraa by default
#    old behaviour avail with -B option
#
#    Revision 1.6  1999/06/02 19:18:21  delwood
#    temp fix for yank and fixed parsing of names with "-"s.
#
#    Revision 1.5  1999/01/12 15:26:54  delwood
#    added ability to read mul format seed files.
#
#    Revision 1.4  1998/11/20 16:59:26  delwood
#    revised search and yank to work on linux or SunOS
#
#    Revision 1.3  1998/09/15 17:15:26  delwood
#    added -p pattern option to search nraa for pattern and add to mini_db
#
#    Revision 1.2  1998/08/26 21:01:56  delwood
#    revised default database, yank_index, and blast options
#    added inspection of existing mini_db to avoid duplicate sequences in db
#    simplified usage
#
#    Revision 1.1  1998/06/03 21:24:07  delwood
#    Initial revision
#
###############################################################

require "getopts.pl";
$usage = <<_EOT_;
Usage:  preseach_db.pl -A fasta_to_search -s mini_db_name
    This searches each member of a fasta file vs db and gets minidatabase of 
    significant hits. It then searches the mini_db with the hmm. 

Options:
 -A fasta	multi-fasta file of seed ("_" or "-" or "." are ignored)
 -R dir         save directory for mini_db, hmm search output 
                    (default is current directory, can use /usr/local/projects/family/tmp)
 -L log_file    write log_file and record search and timing information (not implemented)
 -p pattern     search nraa for pattern (default case indepent) and add these accessions to the mini_db 
                    (it will tell you the names it added this way)
 -i             make nraa search case dependent (ie turn OFF -i flag in grep command)
 -D		debug on
 -h 		help

                New behavior is to search /usr/local/db/panda/nraa/nraa, then yank using the yank server 
 -B             old behavior, uses yank.old, blast on DB (see -d, -Y)
 -d db		database to blast (default /usr/local/projects/family/build_dir/db/nraa_clean)
                (database requires running "setdb" for blasting and "index_to_hash_fasta2" for yanking)
 -Y yank_index  specify if different db (default is "/usr/local/projects/family/build_dir/db/yank_index")


 ____blast options____
 -c cutoff	wu_blastp cutoff e-val for inclusion in db. Default is 10, Maximum is 100
 -o blast_opts  wu_blastp options (default is "-matrix BLOSUM62 -B700 -V700")
                (B and V control number of HSPs and Descriptions returned)
                (add P=1 or P=2 to limit number of processors used)
                (add "-filter "seg+xnu"" to mask low complexity regions - NOT RECOMMENDED)
 -s name	save minidatabase as name for further rounds of searching....
 ____HMM options____
 -H xxx.HMM	name of hmm to be searched
 -P hmmsw	name of hmm searching program
 -O "options" 	hmm searching options (default is "-F -t 10")
 -S hits_name   name to save hmm hits file (default is "hits_hmmname_mini")

added code to check machine type and program location before running. 
Delwood Richardson  Mon Jun  1 14:40:40 1998
based on a script by Jeremy Peterson
_EOT_


&Getopts('H:P:S:O:A:d:Y:R:L:p:iDhs:c:o:xB');
$DEBUG = $opt_D;
die "$usage\n" if $opt_h == 1;
if(length($opt_H) > 0){
    $HMM = $opt_H;
}else{
    print STDERR "no HMM specified, so not searching\n";
    $HMM = "";
}
if(length($opt_A) > 0){
    $SEED = $opt_A;
}else{
    print STDERR "SEED not specified, NO BLAST SEARCH\n";
}
if(length($opt_P) > 0){
    $SEARCH_PROG = $opt_P;
}else{
    print STDERR "HMM Search program not specified, so not searching\n";
}
if(length($opt_p) > 0){
    $PATTERN = $opt_p;
    print STDERR "Searching nraa for pattern \"$PATTERN\", will add to mini_db.\n";
}
if($opt_i == 0){
    $CASE = "-i ";
}else{
    $CASE = "";
}

if(length($opt_O) > 0){
    $SEARCH_OPTS = $opt_O;
}else{
    $SEARCH_OPTS = "-F -t 10";
}
if(length($opt_S) > 0){
    $H_SAVE = $opt_S;
}else{
    $H_SAVE = "";
}
if(length($opt_s) > 0){
    $SAVE = $opt_s;
}else{
    $SAVE = "";
}
if(length($opt_R) > 0){
    $SAVE_DIR = $opt_R;
    $SAVE_DIR =~ s/\/$//;    
}else{
##    $SAVE_DIR = "/usr/local/projects/family/tmp";
    $SAVE_DIR = ".";
}

if(length($opt_c) > 0){
    $E_VAL = $opt_c;
}else{
    $E_VAL = "10";
}
if(length($opt_o) > 0){
    $BLAST_OPTS = $opt_o;
}else{
    $BLAST_OPTS = "-matrix BLOSUM62 -B700 -V700";
}

if($opt_B == 1){
    $OLD_BEHAVE = 1;
    print STDERR "OLD BEHAVIOUR: old yank, " if $DEBUG;
    if(length($opt_d) > 0){
	$DB = $opt_d;
    }else{
	$DB = "/usr/local/projects/family/build_dir/db/nraa_clean";
    }
    if(length($opt_Y) > 0){
	$YANK_INDEX = $opt_Y;
    }else{
	$YANK_INDEX = "/usr/local/projects/family/build_dir/db/yank_index";
    }
}else{
    $DB = "/usr/local/db/panda/nraa/nraa";
    $OLD_BEHAVE = 0;
    print STDERR "NEW BEHAVIOUR: using yank -t, " if $DEBUG;
}
print STDERR "USING DB $DB for search\n" if $DEBUG;



### end options
print STDERR "DEBUG on\n" if $DEBUG;

select STDOUT;
$| = 1;

## $ENV{"BLASTFILTER"} = "/usr/local/util/";
## $ENV{"BLASTMAT"} = "/usr/local/src/ncbi/blastapp./matrix";
$TEMP_DIR = "/usr/local/projects/family/tmp";

##  check if mini_db can be opened (for read/write), otherwise bail!
##  read all current accessions into list

$db_name = ($SAVE eq "" ? "mini_db_$$.fa"  : $SAVE);
$mini_db = "$SAVE_DIR/$db_name";
if(-e $mini_db && -w $mini_db){
    print STDERR "Looking at the current accessions in $mini_db" if $DEBUG;
    ## get the current list out of mini_db
    $list = &check_db($mini_db);
    print STDERR "....done, found $list\n" if $DEBUG;
    $old_len = length($list);
    $old_db_size = scalar(split(/\s+/, $list));
    print STDERR "..... found $old_db_size acc in $mini_db, now searching for more\n";
}elsif(open(TEST, ">$mini_db")){
    close(TEST);
}else{
    die "can't write to $mini_db, so bailing\n$usage\n";
}

## open SEED, make hash of sequences, search db, collect names, yank
## deal with name/coods correctly

$stamp1 = scalar localtime;
goto GREP if($SEED eq "");

open(F_IN, "$SEED") || die "can't open $SEED for reading.\n";
while(<F_IN>){
    chomp;
    if($_ =~ s/^>//){
	if($FASTA == 1){
	    $seqs{$id} = $buf;
	    $id = $buf = "";
	}else{
	    $FASTA = 1;
	}
	($id, $dummy) = split(/\s+/, $_, 2);
	$id =~ s/\/.+//;
    }elsif($FASTA == 1){
	$_ =~ s/[\._-]//g;
	$buf .= $_;
    }else{ ## will take MUL as well. 
	($id, $buf) = split(/\s+/, $_, 2);
	$id =~ s/\/.+//;
	$buf =~ s/[\._-]//g;
	if(length($buf) > 0 && $seqs{$id} eq ""){
	    $seqs{$id} = $buf;
	    $id = $buf = "";
	}else{
	    print STDERR "problem with reading sequence from $id\n";
	}
    }
}
if($FASTA == 1){
    $seqs{$id} = $buf;
    $id = $buf = "";
}
close(F_IN);

@names = keys %seqs;
print STDERR "searching with these sequences: @names\n" if $DEBUG;
foreach $name (@names){
    print STDERR "$name (" . length($seqs{$name}) . " aa)  \t";
    $list = &blast_em($name, $seqs{$name}, $E_VAL, $list, $BLAST_OPTS, $TEMP_DIR);

}

GREP: if($PATTERN ne ""){
    ## search nraa for this pattern and add to list
    print STDERR "Searching for ..$PATTERN.. now\n" if $DEBUG;
    $list = &grep_nraa($list, $PATTERN);
}

$hmm_name = $HMM;
if($hmm_name =~ /\/?([a-zA-Z0-9\-\_]+)\.HMM$/){
    $hmm_name = $1;
}else{
    $hmm_name =~ s/.HMM$//;
}

$cut = length($list) - $old_len + 1;
$new_list = substr($list, $old_len, $cut);
print STDERR "list:..$list..\nnew:..$new_list..\n\n" if $DEBUG;

$new_db_size = &yank_em($new_list, $DB, $YANK_INDEX, $mini_db);
$db_size = $new_db_size + $old_db_size;
print STDERR "$hmm_name FINAL LIST($db_size, added $new_db_size): $new_list\n\n" if $DEBUG;
print STDERR "mini-db of $db_size sequences (added $new_db_size) is at $mini_db\n";

$stamp2 = scalar localtime;

if(length($HMM) > 0 && length($SEARCH_PROG) > 0){
    if($H_SAVE eq ""){
	$output = "$SAVE_DIR/hits_$hmm_name" . "_mini";
    }else{
	$output = "$SAVE_DIR/$H_SAVE";
    }
    &hmm_search($HMM, $SEARCH_PROG, $SEARCH_OPTS, $mini_db, $output);
    $stamp3 = scalar localtime;
    print STDERR "BEGIN: $stamp1, HMM $hmm_name start: $stamp2\nDONE: $stamp3\n";
    print STDERR "search output saved as $output\n" if $DEBUG;
}else{
    print STDERR "BEGIN: $stamp1, DB COMPLETE: $stamp2\n";
}

if($SAVE eq "" && !$DEBUG){
    unlink $mini_db;
}

exit(0);

####Subs
####
sub check_db{
    local($file) = @_;
    local($list, $acc, $dummy);

    open(F_IN, "$file") || die "can't open MINI_DB $file\n";
    while(<F_IN>){
	chomp;
	if($_ =~ s/^>//){
	    ($acc, $dummy) = split(/\s+/, $_, 2);
    	    $list .= "$acc ";
	}
    }
    close(F_IN);
    return($list);
}

sub blast_em{
    local($name, $seq, $E_VAL, $list, $BLAST_OPTS, $TEMP_DIR) = @_;
    local($temp_file, $cmd, @data);
    local($BTAB, $BLAST);
    if($ENV{'HOSTTYPE'} =~ /linux/){
	$BTAB = "/usr/local/bin/btab";
	$BLAST = "/usr/local/bin/blastp";
    }elsif($ENV{'HOSTTYPE'} eq "sun4" && $ENV{'BLASTDIR'} ne ""){
	$BTAB = "/usr/local/util/btab";
	$BLAST = "$ENV{'BLASTDIR'}/blastp";
    }else{
	$BTAB = "/usr/local/util/btab";
	$BLAST = "/usr/local/packages/wu-blast.sparc/blastp     ";
    }
    print STDERR "HOSTTYPE: $ENV{'HOSTTYPE'}\nBLASTDIR: $ENV{'BLASTDIR'}\nBLAST: $BLAST\n" if $DEBUG;

    $temp_file = "$TEMP_DIR/$$.seq";
    open(F_OUT, ">$temp_file") || die "blast_em: can't write seq to $temp_file\n";
    print F_OUT ">$name\n";
    &write_seq(F_OUT, $seq);
    close(F_OUT);
    if($DEBUG){
	print STDERR "searching $name\n$seq\n\n";
    }

    if($E_VAL > 10){
	if($E_VAL <= 100){
	    $E_OPT = "E=$E_VAL ";
	}else{
	    $E_OPT = "E=100 ";
	}
    }else{
	$E_OPT = "";
    }
    $cmd = "$BLAST $DB $temp_file $E_OPT $BLAST_OPTS > $temp_file.nr";
    print STDERR "BLAST: $cmd\n" if $DEBUG;
    system($cmd);
    $cmd = "$BTAB -q $temp_file.nr";
    system($cmd);
    unlink $temp_file if !$DEBUG;
    unlink "$temp_file.nr" if !$DEBUG;
    open(F_IN, "$temp_file.nr.btab") || die "can't open btab file $temp_file.nr.btab";
    while(<F_IN>){
	chomp;
	@data = split /\t/;
	## print STDERR "hit: $data[5]\te_val: $data[19]\tp_val: $data[20]\n" if $DEBUG;
	if($data[19] <= $E_VAL && $list !~ /\Q$data[5]\E /){ # \Q..\E makes | in data[5] literals
	    $list .= "$data[5] ";
	    print STDERR "keeping $data[5] $data[19] $data[20]\n" if $DEBUG;
	}elsif ($list =~ /$data[5] /){
	    print STDERR "skipping $data[5] $data[19] $data[20], already seen\n" if $DEBUG;
	}
    }
    close(F_IN);
    unlink "$temp_file.nr.btab" if !$DEBUG;
    return($list);
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

sub yank_em{
    local($list, $DB, $YANK_INDEX, $mini_db) = @_;
    local(@names, $cmd, $l, $db_size);
    local($yank_dir, $yank_index);

    @names = split(/\s+/, $list);
    $db_size = @names;

    if($OLD_BEHAVE == 1){
	$yank_dir = $YANK_INDEX;
	$yank_dir =~ s/(\w+)$//;
	$yank_index = $1;
	if($ENV{'HOSTTYPE'} =~ /linux/){
	    $yank_index .= "_linux";
	}


	for($i=0;$i < $db_size; $i+= 5){
	    $l = "$names[$i] $names[$i+1] $names[$i+2] $names[$i+3] $names[$i+4] ";
	    $cmd = "echo '$l' | /usr/local/bin/yank.old -I $yank_dir -i $yank_index  >> $mini_db";
##	$cmd = "echo '$l' | yank -Y -I $yank_dir -i $yank_index  >> $mini_db";
	    print STDERR "yank cmd: $cmd\n" if $DEBUG;
	    system $cmd;
	}
    }else{
	for($i=0;$i < $db_size; $i+= 5){
	    $l = "$names[$i] $names[$i+1] $names[$i+2] $names[$i+3] $names[$i+4] ";
	    $cmd = "echo '$l' | /usr/local/bin/yank -t  >> $mini_db";
	    print STDERR "yank cmd: $cmd\n" if $DEBUG;
	    system $cmd;
	}
    }
    return($db_size);
}

sub hmm_search{
    local($HMM, $SEARCH_PROG, $SEARCH_OPTS, $mini_db, $output) = @_;
    local($cmd);
    
    $cmd = "$SEARCH_PROG $SEARCH_OPTS $HMM $mini_db > $output";
    print STDERR "HMM search: $cmd\n" if $DEBUG;
    system $cmd;

}

sub grep_nraa{
    local($list, $PATTERN) = @_;
    local($new_list, @return, $acc, $def, $dummy, $cmd);
    local($NRAA) = "/usr/local/db/panda/nraa/nraa";

    $new_list = $list;
    $cmd = "grep \">\" $NRAA | grep $CASE \"$PATTERN\" ";
    print STDERR "GREP CMD: $cmd\n" if $DEBUG;
    open(C_IN, "$cmd |") || die "can't seem to grep nraa using '$cmd'\n";
    while(<C_IN>){
	chomp;
	($def, $dummy) = split(/\cA/, $_, 2);
	if($def ne ""){
	    ($acc, $dummy) = split(/\s+/,$def, 2);
	    $acc =~ s/^>//;
	    print STDERR "acc:..$acc.. DEF:..$def..\n" if $DEBUG;
	    if($list !~ /\Q$acc\E /){
		print STDERR "NEW: $def\n";
		$new_list .= "$acc ";
	    }else{
		print STDERR "IN DB: $def\n";
	    }
	}
    }
    close(C_IN);
    return($new_list);
}

 
