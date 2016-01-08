#!/usr/bin/env perl -w

# Author: Fayaz Seifuddin, Nora Khaldi
# Date Created: 06/05/2007
# Last Date Modified: 06/05/2007
# File:testscript4.1.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: Read in CLUSTERINFOOUT and HMMoutput.htab (parsed from htab.pl) files to compare and get necessary information, STEP4
# Domain-content-clusters
# Usage: perl testscript4.1.pl GENEDISTOUT $tempdir

use strict;

# Varible declaration
my $filename;  # Name of the input file1
my $filename1; # Names of domains/files from HMMoutput.htab file (from htab.pl)
my %hash1;     # Hash1 to store header and sequence information from CLUSTERINFOOUT file
my $counter=-1; # Counter to loop through the CLUSTERINFOOUT file
my $line; # Each line from CLUSTERINFOOUT file
my @headers; # Array to store headers from CLUSTERINFOOUT file
my @sequences; # Array to store sequences from CLUSTERINFOOUT file
my $i; # Counter to create hash from the CLUSTERINFOOUT file
my $k; # store enzyme_id (NRPS,PKS,DMAT ids)for hash2
my $l; # store gene_id for hash2
my $m; # store position of gene with respect to core enzyme (+/- 20) 
my $tempdir; # folder in which all HMMsearch outputs are stored
#my $b=-1; # Counter to iterate through all domains/files from HMMoutput2.htab file (from htab.pl)

# Read from CLUSTERINFOOUT file containing multiple sequences and headers (gene clusters)
# Open the files given as arguments on the command line
$filename = $ARGV[0];  # inputFASTA file
$tempdir = $ARGV[1]; # $tempdir folder in which all HMMseach outputs are stored

# write cluster-domain-content-out to output file
#open OUTFILE,"> CLUSTER-DOMAIN-CONTENT-OUT" or die "Can't open: OUT\n";

open(INFILE, $filename) or die "Can't open $filename";
# Read in headers and sequencs from file1
while ($line = <INFILE>){
    chomp $line;                         # Remove extra end line \n from EOF
    # If loop to skip to the next line if the first line in file is blank
    if ($line =~ /^\s*$/){
        next;
    }
    elsif ($line =~ /^>/) {               # Line starts with a ">"
        $counter++;                       # Increament counter to go to next line
        $headers[$counter]=($line);       # Save header line
        $sequences[$counter] = '';        # Start with an empty sequence
    }
    else {
        $sequences[$counter]=$sequences[$counter].$line; # Concatenate sequence
    }
}# End while loop
close(INFILE);

# Loop to create hash1 from CLUSTERINFOOUT file headers
for (my $i=0; $i<=scalar(@headers)-1; $i++){
    my @newelement = split(/\t/,$headers[$i]); # split headers by tab
    my $enzyme_id = $newelement[0]; # store enzyme_id
    my $gene_id = $newelement[1]; # store gene_id
    my $position = $newelement[2]; # store position of gene from core enzyme
    #ABC_membrane, ABC_tran.out, Acetyltransf_1, adh_short, ADH_zinc_N,FAD_binding_3,FAD_binding_4,
    #Fungal_trans,hmm_AA_permease,Methyltransf_2,MFS_1,p450,Pyr_redox,Sugar_tr,Transferase,Zn_clus
    $hash1{$enzyme_id}{$gene_id}{$position} = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]; # add domain content info as arrays
    $hash1{$enzyme_id}{$gene_id}{$position}[31]  = $newelement[4]; # gene_function
    $hash1{$enzyme_id}{$gene_id}{$position}[32]  = $newelement[5]; # contig/chromosome
    $hash1{$enzyme_id}{$gene_id}{$position}[33]  = $newelement[6]; # gene_order
    $hash1{$enzyme_id}{$gene_id}{$position}[34]  = $newelement[7]; # gene_start5
    $hash1{$enzyme_id}{$gene_id}{$position}[35]  = $newelement[8]; # gene_end3
    $hash1{$enzyme_id}{$gene_id}{$position}[36]  = $newelement[3]; # gene_distance
    $hash1{$enzyme_id}{$gene_id}{$position}[37]  = $sequences[$i]; # protein sequences
}# End for loop

# Read in .htab files from /HMMoutput2 directory
my @files = <$tempdir/*.output.htab>;

# Loop through each HMMoutput2.htab file and get the necessary information
foreach $filename1(@files){
    #$b=$b+1; # increment counter to iterate to next file
    my @r = []; # store information from each file (each column)
    my @r1 = []; # store information from biological function, get gene_id & position
    my $enzyme_id = ''; # store enzyme_id
    my %hash2; # store domain content in hash
    #print OUTFILE "$filename1\n"; # to determine the order of the domains i.e. order of columns in table.
    #print "$filename1\n";
    # Read in HMMouput file (from htab.pl), create hash2
    local $/="\n";  # going back to default setting '\n' being the seperator
    open(INFILE, $filename1) or die "Can't open $filename1";
    while (<INFILE>) {
        @r=split(/\t/,$_); # split the file by \t
        chomp @r; # remove extra \n
        $enzyme_id = ">".$r[5]; # concatenate ">" infront of gene_id to match
	@r1=split(/ /,$r[16]); # split biological function by ' ',get gene_id & position of gene with respect to core enzyme
        $hash2{$enzyme_id}{$r1[0]}{$r1[1]}[0] = $r[14]; # total number of domains as read from htab.pl
    }#end while

    # Loop to check if gene_ids are the same in both files (inputFASTA file & HMMoutput)
    # & if the same add the info. from the first file to that
    # gene_id--->to prevent mismatches
    # create hash1 (hash containing all the information)
    foreach $k(keys %hash1){ # enzyme_id
	foreach my $l(keys %{$hash1{$k}}){ # gene_id
	    foreach my $m(keys %{$hash1{$k}{$l}}){ # position of gene with respect to core enzyme
		if(exists $hash2{$k}{$l}{$m}){ # check to see if it matches with CLUSTERINFOOUT to prevent mismatches
		    #$hash1{$k}{$l}{$m}[$b] = $hash2{$k}{$l}{$m}[0]; # add total number of domains as an array
		    if($filename1 eq "$tempdir/2OG-FeII_Oxy.output.htab"){
			$hash1{$k}{$l}{$m}[0] = $hash2{$k}{$l}{$m}[0];
		    }
		    if($filename1 eq "$tempdir/ABC_membrane.output.htab"){
                        $hash1{$k}{$l}{$m}[1] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/ABC_tran.output.htab"){
                        $hash1{$k}{$l}{$m}[2] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Acetyltransf_1.output.htab"){
                        $hash1{$k}{$l}{$m}[3] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/adh_short.output.htab"){
                        $hash1{$k}{$l}{$m}[4] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/ADH_zinc_N.output.htab"){
                        $hash1{$k}{$l}{$m}[5] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Aminotran_1_2.output.htab"){
                        $hash1{$k}{$l}{$m}[6] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/FAD_binding_3.output.htab"){
                        $hash1{$k}{$l}{$m}[7] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/FAD_binding_4.output.htab"){
                        $hash1{$k}{$l}{$m}[8] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Fungal_trans.output.htab"){
                        $hash1{$k}{$l}{$m}[9] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/GST_C.output.htab"){
                        $hash1{$k}{$l}{$m}[10] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/GST_N.output.htab"){
                        $hash1{$k}{$l}{$m}[11] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/GTP_cyclohydroI.output.htab"){
                        $hash1{$k}{$l}{$m}[12] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/hmm_AA_permease.output.htab"){
                        $hash1{$k}{$l}{$m}[13] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Lactamase_B.output.htab"){
                        $hash1{$k}{$l}{$m}[14] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/methyl_EasF.output.htab"){
                        $hash1{$k}{$l}{$m}[15] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Methyltransf_2.output.htab"){
                        $hash1{$k}{$l}{$m}[16] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/MFS_1.output.htab"){
                        $hash1{$k}{$l}{$m}[17] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/p450.output.htab"){
                        $hash1{$k}{$l}{$m}[18] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Peptidase_M19.output.htab"){
                        $hash1{$k}{$l}{$m}[19] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Peptidase_M24.output.htab"){
                        $hash1{$k}{$l}{$m}[20] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/PhyH.output.htab"){
                        $hash1{$k}{$l}{$m}[21] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Pyr_redox.output.htab"){
                        $hash1{$k}{$l}{$m}[22] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Pyr_redox_2.output.htab"){
                        $hash1{$k}{$l}{$m}[23] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Sugar_tr.output.htab"){
                        $hash1{$k}{$l}{$m}[24] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Transferase.output.htab"){
                        $hash1{$k}{$l}{$m}[25] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Zn_clus.output.htab"){
                        $hash1{$k}{$l}{$m}[26] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Epimerase.output.htab"){
                        $hash1{$k}{$l}{$m}[27] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/ergot_EASG.output.htab"){
                        $hash1{$k}{$l}{$m}[28] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/Oxidored_FMN.output.htab"){
                        $hash1{$k}{$l}{$m}[29] = $hash2{$k}{$l}{$m}[0];
                    }
		    if($filename1 eq "$tempdir/tRNA-synt_2.output.htab"){
                        $hash1{$k}{$l}{$m}[30] = $hash2{$k}{$l}{$m}[0];
                    }
		}#end if loop
	    }#end foreach loop-{%hash1}
	}#end foreach loop-{$hash1{$k}}
    }#end foreach loop-{$hash1{$k}{$l}}
close(INFILE);
}# end foreach

# write cluster-domain-content-out to output file
#open OUTFILE,"> CLUSTER-DOMAIN-CONTENT-OUT" or die "Can't open: OUT\n";

# print hash2, cluster information & domain content information
foreach my $enzyme(sort keys %hash1){ # enzyme_id
    foreach my $gid(sort keys %{$hash1{$enzyme}}){ # gene_id
	foreach my $loc(sort keys %{$hash1{$enzyme}{$gid}}){ # position of gene with respect to core enzyme
	    #print "$enzyme\t$gid\t$loc\t$hash1{$enzyme}{$gid}{$loc}[16]\t$hash1{$enzyme}{$gid}{$loc}[17]\t$hash1{$enzyme}{$gid}{$loc}[18]\t$hash1{$enzyme}{$gid}{$loc}[19]\t$hash1{$enzyme}{$gid}{$loc}[20]\n$hash1{$enzyme}{$gid}{$loc}[21]\n";
	    #print "$enzyme\t$gid\t$loc\n";

	    #if($hash1{$enzyme}{$gid}{$loc}[0]==0 && $hash1{$enzyme}{$gid}{$loc}[1]==0 && $hash1{$enzyme}{$gid}{$loc}[2]==0 && $hash1{$enzyme#}{$gid}{$loc}[3]==0 && $hash1{$enzyme}{$gid}{$loc}[4]==0 && $hash1{$enzyme}{$gid}{$loc}[5]==0 && $hash1{$enzyme}{$gid}{$loc}[6]==0 && $hash1{#$enzyme}{$gid}{$loc}[7]==0 && $hash1{$enzyme}{$gid}{$loc}[8]==0 && $hash1{$enzyme}{$gid}{$loc}[9]==0 && $hash1{$enzyme}{$gid}{$loc}[10]==0 &&#$hash1{$enzyme}{$gid}{$loc}[11]==0 && $hash1{$enzyme}{$gid}{$loc}[12]==0 && $hash1{$enzyme}{$gid}{$loc}[13]==0 && $hash1{$enzyme}{$gid}{$loc}#[14]==0 && $hash1{$enzyme}{$gid}{$loc}[15]==0 && $hash1{$enzyme}{$gid}{$loc}[16]==0 && $hash1{$enzyme}{$gid}{$loc}[17]==0 && $hash1{$enzyme}{$gid}{$loc}[18]==0 && $hash1{$enzyme}{$gid}{$loc}[19]==0 && $hash1{$enzyme}{$gid}{$loc}[20]==0 && $hash1{$enzyme}{$gid}{$loc}[21]==0 && $hash1{$enzyme}{$gid}{$loc}[22]==0 && $hash1{$enzyme}{$gid}{$loc}[23]==0 && $hash1{$enzyme}{$gid}{$loc}[24]==0 && $hash1{$enzyme}{$gid}{$loc}[25]==0 && $hash1{$enzyme}{$gid}{$loc}[26]==0){ # eliminate all genes without a domain hit
	    #}#end if
	    #else{
	    print "$enzyme\t$gid\t$loc\t$hash1{$enzyme}{$gid}{$loc}[32]\t$hash1{$enzyme}{$gid}{$loc}[33]\t$hash1{$enzyme}{$gid}{$loc}[34]\t$hash1{$enzyme}{$gid}{$loc}[35]\t$hash1{$enzyme}{$gid}{$loc}[36]\t$hash1{$enzyme}{$gid}{$loc}[0]\t$hash1{$enzyme}{$gid}{$loc}[1]\t$hash1{$enzyme}{$gid}{$loc}[2]\t$hash1{$enzyme}{$gid}{$loc}[3]\t$hash1{$enzyme}{$gid}{$loc}[4]\t$hash1{$enzyme}{$gid}{$loc}[5]\t$hash1{$enzyme}{$gid}{$loc}[6]\t$hash1{$enzyme}{$gid}{$loc}[7]\t$hash1{$enzyme}{$gid}{$loc}[8]\t$hash1{$enzyme}{$gid}{$loc}[9]\t$hash1{$enzyme}{$gid}{$loc}[10]\t$hash1{$enzyme}{$gid}{$loc}[11]\t$hash1{$enzyme}{$gid}{$loc}[12]\t$hash1{$enzyme}{$gid}{$loc}[13]\t$hash1{$enzyme}{$gid}{$loc}[14]\t$hash1{$enzyme}{$gid}{$loc}[15]\t$hash1{$enzyme}{$gid}{$loc}[16]\t$hash1{$enzyme}{$gid}{$loc}[17]\t$hash1{$enzyme}{$gid}{$loc}[18]\t$hash1{$enzyme}{$gid}{$loc}[19]\t$hash1{$enzyme}{$gid}{$loc}[20]\t$hash1{$enzyme}{$gid}{$loc}[21]\t$hash1{$enzyme}{$gid}{$loc}[22]\t$hash1{$enzyme}{$gid}{$loc}[23]\t$hash1{$enzyme}{$gid}{$loc}[24]\t$hash1{$enzyme}{$gid}{$loc}[25]\t$hash1{$enzyme}{$gid}{$loc}[26]\t$hash1{$enzyme}{$gid}{$loc}[27]\t$hash1{$enzyme}{$gid}{$loc}[28]\t$hash1{$enzyme}{$gid}{$loc}[29]\t$hash1{$enzyme}{$gid}{$loc}[30]\t$hash1{$enzyme}{$gid}{$loc}[31]\n";
	    #}#end else   
	}#end foreach {$hash1{$enzyme}{$gid}}
    }#end foreach {$hash1{$enzyme}}
}#end foreach {%hash1}

#close OUTFILE;

exit;










