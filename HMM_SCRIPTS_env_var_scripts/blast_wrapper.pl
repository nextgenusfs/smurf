#!/usr/local/bin/perl
## questions: remove "redundant" sequences that have identical accession??
{
    use strict;
    use Cwd;
    require "$ENV{HMM_SCRIPTS}/hmm_lib.pl"; #???
#    require "/usr/local/projects/OGC/egad_test/working/hmm_lib.pl";

    my $usage ="
NAME
    blast_wrapper.pl - run wu-blast(or NCBI blastall??, etc) against a database and
                       accumulate results obtained.

SYNOPSIS
    blast_wrapper.pl [options]

DESCRIPTION
 run blast:
    -B   blast_command (default: /usr/local/common/blastp)
         string for running blast

    -D   db_name (default: /usr/local/db/panda/nraa/nraa)
         database name against that blast searches

    -O   blast_options (default: -matrix BLOSUM62 -B1000 -V1000)
         string of blast option

    -S   sequence_file (default: null)
         sequence file including its path for expanding mini_db or changing accession

 job control
    -e   expand|new for mini_db_mode (default: no calcualtion)
         expand : create a new or expand    an old mini_db file specified by -t
         new:     create a new or overwrite an old mini_db file specified by -t

    -a   reorder|convert (default: no change to protein accession)
         reorder: rename protein accession in the order EGAD -> OMNI -> SP -> GP -> PIR
                  if the original protein accession can be found by YANK server. No
                  request for blastp in this mode.
         convert: convert the original protein accession to an entry in nraa that has
                  the inout sequence as its substring or with high similarity.  May request
                  for blastp if necessary.

 general setup
    -d   directory (default: current directory)
         directory where all results such as new or expanded mini_db and/or new sequence file
         are saved.  It is also used for saving temporary files.

    -t   mini_db_name (default: ./mini_db)
         name of mini_db file including its path. Note that new or expanded mini_db file is
         saved in the directory specified by -d option.  The original mini_db file is not
         modified if it is not under the directory specified by -d option

    -o   new_seq_file (default: new_seq_file)
         name of output sequence file with the original accession being replaced by the first
         accession in the database.  This file is saved in the directory specified by -d option. 

    -f   format_output (default: fasta)
         output sequence file format for new_seq_file. Valid options are: fasta, mul, msf.

    -l   logfile_name (default: sequence_file name followed by .blast_wrapper_log)
         It is save in the directory specified by -d option.

    -h   print this message";
    
    my %info;
    my $patt1 = '^([^\|\:\.\/\s]+((\||\:|\.)[^\|\/\s\.\:]+)?)';
    &get_options(\%info, $usage);
    unlink "$info{log_file}";
    #prepare mini_db file in save_dir
    my $mini_db_file_name = $1 if($info{mini_db_file} =~ /([^\/\s]+)$/);
    $info{new_mini_db_file} = "$info{save_dir}/$mini_db_file_name";
    if($info{mini_db_mode} eq "expand") {
	system "cp $info{mini_db_file} $info{new_mini_db_file}";
    }
    elsif($info{mini_db_mode} eq "new") {
	unlink "$info{new_mini_db_file}";
    }
    #read seq file and/or mini_db file
    my (%seq, %mini_seed, %new_seq);
    if(-s $info{seq_file}) {
	system "chmod 444 $info{seq_file}" unless(-r $info{seq_file});
	&read_sequence(\%seq, $info{seq_file}, '', '', '') if(-s $info{seq_file});
	%new_seq = %seq;
    }
    if(-s $info{new_mini_db_file}) {
	system "chmod 444 $info{new_mini_db_file}" unless(-r $info{new_mini_db_file});
	&read_sequence(\%mini_seed, $info{new_mini_db_file}, '', '', '') if(-s $info{new_mini_db_file});
    }
    #get the short acc for each prot in seq
    my %tmp_yank;
    for my $n (0 .. $seq{number} - 1) {
#	if($seq{$n}{header} =~ /$patt1/){
	if($seq{$n}{header} =~ /^([^\|]+\|[^\|\/\s]+)/){
	    $tmp_yank{$1} = {};
	    $seq{$n}{short_acc} = $1; #this is needed for general purpose
	}
	else {
	    &print_cyc($info{log_file}, "${\(scalar localtime)}: ***Error: unusual protein accession($seq{$n}{header}) causes incomplete construction of mini_db\n", 0, 1);
	}
    }
    &yank_prot2(\%tmp_yank, '', '') if($info{convert_acc} ne ''); # yank only if protein accessions need modification
    #main process here
    for my $n (0 .. $seq{number} - 1) {
	my ($run_blast, $in_seed, $converted) = (0, 0, 0);
	$converted = 1 if($seq{$n}{header} =~ /^SP/ or $info{convert_acc} eq "");
	#determin whether the sequence is part of the seed of mini_db(blastp is unnecessary if yes)
	if($info{mini_db_mode} ne "") {
	    for my $m (0 .. $mini_seed{stockholm}{GF}{mini_db}{number} - 1) {
		$in_seed = 1 if($mini_seed{stockholm}{GF}{mini_db}{$m}{seq} =~ /\Q$seq{$n}{seq}\E/i);
	    }
	}
	#reorder protein accession
	if($info{convert_acc} =~ /^(reorder|convert)$/ && $converted == 0 && $seq{$n}{short_acc} ne '' && $tmp_yank{$seq{$n}{short_acc}}{seq} ne '' && $tmp_yank{$seq{$n}{short_acc}}{header} ne '') {
	    if($tmp_yank{$seq{$n}{short_acc}}{seq} =~ /\Q$seq{$n}{seq}\E/i) {
		if($tmp_yank{$seq{$n}{short_acc}}{header} =~ /^\s*(\S+)/) {
		    my $tmp = $1." ".$seq{$n}{seq_gap};
		    $new_seq{$n} = {};
		    $new_seq{$n}{ori} = $tmp;
		    $new_seq{$n}{header} = $1;
		    &parse_single_seq($new_seq{$n}, 'mul', '', '', '', '', '');
		    $converted = 1;
		}
		else {
		    &print_cyc($info{log_file}, "${\(scalar localtime)}: +++Error: empty header for $seq{$n}{short_acc}\n", 0, 1);
		}
	    }
	    else {
		&print_cyc($info{log_file}, "${\(scalar localtime)}: +++Error: the sequence (\#$n) retrieved from DB is not identical to the original one\n", 0, 1);
	    }
	}
	#convert protein accession if its sequence is seen before
	if($info{convert_acc} eq "convert" and $converted == 0) {
	    for my $m (0 .. $mini_seed{number}-1) {
		my $match_seq = $seq{$n}{seq};
		#further process $match_seq eg remove initial Met etc.
		if($mini_seed{$m}{seq} =~ /\Q$match_seq\E/i) { # this condition can be replaced by blastp match and identity > eg. 98%
		    my $tmp = "$mini_seed{$m}{header} $seq{$n}{seq_gap}";
		    $new_seq{$n} = {};
		    $new_seq{$n}{ori} = $tmp;
		    $new_seq{$n}{header} = "$mini_seed{$m}{header}";
		    &parse_single_seq($new_seq{$n}, 'mul', '', '', '', '', '');
		    $converted = 1;
		}
	    }
	}
	# determine whether to run blast or not
	$run_blast = 1 if($in_seed == 0 && $info{mini_db_mode} ne '');
	$run_blast = 1 if($converted == 0  && $info{convert_acc} eq "convert");
	#run blast
	if($run_blast == 1) {
	    my $out = &format_sequence(\%seq, 'fasta', $n, 60, 0, '', '', '');
	    my $tmp_acc = $seq{$n}{short_acc};
	    $tmp_acc =~ s/[^\w]/_/g;
	    my $tmp_file = &tmp_file($info{save_dir}, "tmp_${tmp_acc}_");
	    open(FH, ">$info{save_dir}/$tmp_file.fa");
	    print FH $out;
	    close FH;
	    my $blast_cmd;
	    $blast_cmd = "$info{blast_cmd} $info{db_file} $info{save_dir}/$tmp_file.fa $info{blast_opt} > $info{save_dir}/$tmp_file.log";
	    my $result_blast = system "$blast_cmd";
	    if($result_blast == 256) {
		exit 102; # 256 X 102 = 26112
	    }
	    elsif($result_blast != 0) {
		print "${\(scalar localtime)}: ***Error: unknown error(exit code $result_blast. code 2 or 15, 512 might be from killing of blastp process)\n";
		die;
	    }
	    #process blast results
	    my ($success, $percent, $header_new, $lend, $rend) = &parse_blast("$info{save_dir}/$tmp_file.log", $info{save_dir});
	    #build mini_db
	    if($info{mini_db_mode} ne '' && $in_seed == 0) {
		my ($result, $log) = &build_mini_db($seq{$n}{seq}, $new_seq{$n}{header}, \%mini_seed, $info{db_file}, "$info{save_dir}/$tmp_file.log.btab", $blast_cmd, $info{save_dir});
		print $log if($result eq 'error');
	    }
	    # convert accession
	    if($info{convert_acc} eq "convert" && $converted == 0 && $success == 1) {
		my ($t1) = split(/\s/, $header_new);
		my $t2 = $t1."\t".$seq{$n}{seq_gap};
		%{$new_seq{$n}} = ();
		$new_seq{$n}{ori} = $t2;
		my $fatal = '';
		&parse_single_seq($new_seq{$n}, 'mul', '', \$fatal, '', '', '');
		print "${\(scalar localtime)}: ***Error: $fatal\n" if($fatal ne '');
		$new_seq{$n}{lend} = $lend;
		$new_seq{$n}{rend} = $rend;
		$converted = 1;
	    }
	    unlink "$info{save_dir}/$tmp_file", "$info{save_dir}/$tmp_file.fa", "$info{save_dir}/$tmp_file.log", "$info{save_dir}/$tmp_file.log.btab";
	}
    }
    #output mini_db file
    open(FH, ">$info{new_mini_db_file}");
    my $out = &format_sequence(\%mini_seed, 'mul', '', 60, 0, '', '', '');
    print FH $out;
    close FH;
    #output converted sequences
    if($info{convert_acc} ne '') {
	&get_coor_from_nraa("/usr/local/db/panda/AllGroup/AllGroup.niaa", \%new_seq);
#	&get_coor_from_nraa("/usr/local/db/panda/nraa/nraa", \%new_seq);
	open(FH, ">$info{save_dir}/$info{new_seq_file}");
	my $out = &format_sequence(\%new_seq, $info{new_seq_format}, '', 60, 1, '', '', '');
	print FH $out;
	close FH;
    }
    exit 101; # 101 X 256 = 25856
}
sub build_mini_db{
# $seq and $header are the ungapped sequence (can be fragment) and accession of the protein used to build or expand mini_db
# $mini_seed_r is the hash reference for mini_db
    use strict;
    my($seq, $header, $mini_seed_r, $db_file, $btab_file, $blast_cmd, $save_dir) = @_;
    my $sub_status = "success";
    my $sub_log;
    my $patt1 = '^([^\|\:\.\/\s]+((\||\:|\.)[^\|\/\s\.\:]+)?)';
    #add the seq to mini_seed seed record.
    $$mini_seed_r{stockholm}{GF}{mini_db}{number} = 0 if($$mini_seed_r{stockholm}{GF}{mini_db}{number} eq '');
    my $num = $$mini_seed_r{stockholm}{GF}{mini_db}{number};
    $$mini_seed_r{stockholm}{GF}{mini_db}{$num}{acc}     = $header;
    $$mini_seed_r{stockholm}{GF}{mini_db}{$num}{seq}     = $seq;
    $$mini_seed_r{stockholm}{GF}{mini_db}{$num}{cmd}     = $blast_cmd;
    $$mini_seed_r{stockholm}{GF}{mini_db}{$num}{db_name} = $db_file;
    $$mini_seed_r{stockholm}{GF}{mini_db}{$num}{db_size} = -s $db_file;
    open(FH, "stat $db_file |");
    while(<FH>) {
	$$mini_seed_r{stockholm}{GF}{mini_db}{$num}{db_date} = $1 if($_ =~ /^Modify: (.+)\(/);
    }
    close FH;
    ++$$mini_seed_r{stockholm}{GF}{mini_db}{number};
    # get btab info
    open(FH, "$btab_file");
    my @hits = <FH>;
    close FH;
    my (%tmp_head, %prot, %mini_seed_header); # %tmp_head is used to remove redundant prot acc due to multiple hits to a prot.
    for my $n (0 .. @hits-1) {
	my @tmp = split(/\t/, $hits[$n]);
	$prot{$1} = {} if($tmp[5] =~ /^([^\|]+\|[^\|\/\s]+)/);
#	$prot{$1} = {} if($tmp[5] =~ /$patt1/); # $tmp[5] is the acc of matches
    }
    # add seed seq to hits
    $prot{$1} = {} if($header =~ /$patt1/);
    &yank_prot2(\%prot, '', '');
    for my $n (0 .. $$mini_seed_r{number}-1) {
	$mini_seed_header{$1} = 1 if($$mini_seed_r{$n}{header} =~ /$patt1/);
    }
    my $acc_exist = 0;
    for my $k (keys %prot) {
	if($prot{$k}{seq} ne '' and $prot{$k}{first_acc} =~ /$patt1/){
	    my $acc = $1;
	    my $found = 0;
	    $found = 1 if($mini_seed_header{$acc} == 1);
	    if($found == 0) {
		my %tmp_seq;
		if ($prot{$k}{seq} =~ /$seq/ ){
		    $acc_exist = 1;
		}
		$tmp_seq{ori} = "$prot{$k}{first_acc} $prot{$k}{seq}";
		&parse_sequence(\%tmp_seq, 'mul', '', '', '', '', '');
		%{$$mini_seed_r{$$mini_seed_r{number}}} = %{$tmp_seq{0}};
		++$$mini_seed_r{number};
		$mini_seed_header{$acc} = 1;
	    }
	}
    }
    #put the seed sequence into mini_db if it is not there yet
#    if (!$acc_exist){

    my $acc = $1 if($header =~ /$patt1/);
    if(($mini_seed_header{$acc} == 0) && (!$acc_exist)) {
	my %tmp_seq;
	$tmp_seq{ori} = "$header $seq";
	&parse_sequence(\%tmp_seq, 'mul', '', '', '', '', '');
	%{$$mini_seed_r{$$mini_seed_r{number}}} = %{$tmp_seq{0}};
	++$$mini_seed_r{number};
    }
    return ($sub_status, $sub_log);
}
sub get_options{
    use strict;
    use Getopt::Std;
    use vars qw($opt_B $opt_D $opt_O $opt_S $opt_a $opt_d $opt_e $opt_f $opt_h $opt_l $opt_o $opt_p $opt_t $opt_u);
    my ($info_r, $usage) = @_;
    getopts('B:D:O:S:a:d:e:f:hl:o:p:t:u:') or die "Wrong input options. \n"; #get options
    die "$usage\n" if($opt_h);
    #get blast command, database for search, options, sequence file --------------------------
    $$info_r{blast_cmd} = ($opt_B ? $opt_B : "/usr/local/common/blastp");
    $$info_r{db_file}   = ($opt_D ? $opt_D : "/usr/local/db/panda/AllGroup/AllGroup.niaa");
#    $$info_r{db_file}   = ($opt_D ? $opt_D : "/usr/local/db/panda/nraa/nraa");
    $$info_r{blast_opt} = ($opt_O ? $opt_O : "-cpus 1 -matrix BLOSUM62 -B1000 -V1000");
    if($opt_S) {
	$$info_r{seq_file} = $opt_S;
	#check existance and permission of the sequence file
	unless(-s $$info_r{seq_file}) {
	    print "${\(scalar localtime)}: ***Error: the sequence file $$info_r{seq_file} contains no data or does not exist\n";
	    exit 102;
	}
	system "chmod 664 $$info_r{seq_file}";
	unless(-r $$info_r{seq_file}) {
	    print "${\(scalar localtime)}: ***Error: no read permission for the sequence file $$info_r{seq_file}\n";
	    exit 102;
	}
    }
    die "++Warning: no sequence file specified\n$usage\n" unless($$info_r{seq_file});
    $$info_r{new_seq_file} = ($opt_o ? $opt_o : "new_seq_file");
    $$info_r{new_seq_format} = ($opt_f ? $opt_f : "fasta");
    #create dir for saving files and check its existance
    $$info_r{save_dir} = ($opt_d ? $opt_d : "./");
    if($$info_r{save_dir} ne './'){
	mkdir $$info_r{save_dir}, 0777;
	system "chmod 775 $$info_r{save_dir}";
	unless(-d $$info_r{save_dir}) {
	    print "${\(scalar localtime)}: ***Error: the directory $$info_r{save_dir} for saving files does not exist or can not be created\n";
	    exit 102;
	}
    }
    unless(-r $$info_r{save_dir} && -w $$info_r{save_dir} && -x $$info_r{save_dir}) {
	print "+++Warning: the directory $$info_r{save_dir} does not have enough permission\n";
	#exit 102;
    }
    #mode for mini_db
    if($opt_e =~ /^(expand|new)$/) {
	$$info_r{mini_db_mode} = $opt_e;
    }
    elsif($opt_e ne '') {
	print "${\(scalar localtime)}: ***Error: unknown mode for mini_db calculation\n";
	exit 102;
    }
    else{
	$$info_r{mini_db_mode} = '';
    }
    if($$info_r{mini_db_mode} ne '') {
	$$info_r{mini_db_file} = ($opt_t ? $opt_t : "mini_db");
	if(-e $$info_r{mini_db_file}) {
	    system "chmod 666 $$info_r{mini_db_file}" unless(-r $$info_r{mini_db_file});
	    unless(-r $$info_r{mini_db_file}) {
		print "${\(scalar localtime)}: ***Error: no permission to read the mini_db file $$info_r{mini_db_file}\n";
		exit 102;
	    }
	}
    }
    #mode for protein accession
    if($opt_a =~ /^(reorder|convert)$/) {
	$$info_r{convert_acc} = $opt_a;
    }
    elsif($opt_a ne '') {
	print "${\(scalar localtime)}: ***Error: unknown mode for converting protein accession\n";
	exit 102;
    }
    else{
	$$info_r{convert_acc} = '';
    }
    #
    my $seq_file = $1 if($$info_r{seq_file} =~ /\/?([^\/]+)$/);
    $$info_r{log_file_name} = ($opt_l ? $opt_l : "$seq_file.blast_wrapper_log");
    $$info_r{log_file} = "$$info_r{save_dir}/$$info_r{log_file_name}";
}
