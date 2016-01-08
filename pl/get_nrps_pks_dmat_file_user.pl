#!/usr/bin/env perl -w

# Author: Fayaz Seifuddin, Nora Khaldi
# Date Created: 05/23/2007
# Last Date Modified: 06/04/2007
# File:testscript3.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: Read in inputFASTA file and HMMoutput.htab (parsed from htab.pl) files to compare and get necessary information, STEP3
# Domain-content, determine NRPS,PKS,NRPS-like,PKS-like,DMAT
# Usage: perl testscript3.pl inputFASTA $tempdir

use strict;

# Varible declaration
my $filename;  # Name of the input file1
my $filename1; # Names of domains/files from HMMoutput.htab file (from htab.pl)
my %hash1;     # Hash1 to store header and sequence information from inputFASTA file
my $counter=-1; # Counter to loop through the inputFASTA file
my $line; # Each line from inputFASTA file
my @headers; # Array to store headers from inputFASTA file
my @sequences; # Array to store sequences from inputFASTA file
my $i; # Counter to create hash from the inputFASTA file
my $k; # Gene_id variable for hash2
my $tempdir; # folder for HMMsearch output for each genome

# Read from a inputFASTA file containing multiple sequences and headers
# Open the files given as arguments on the command line
$filename = $ARGV[0];  # inputFASTA file
$tempdir = $ARGV[1]; # folder for HMMsearch output for each genome

# Write domain content out to file
#open OUTFILE,"> DOMAINCONTENTOUT.fasta" or die "Can't open: OUT\n";

open(INFILE, $filename) or die "Can't open $filename";
# Read in headers and sequences from file1
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

# Loop to create hash1 from inputFASTA file headers
for (my $i=0; $i<=scalar(@headers)-1; $i++){
    my @newelement = split(/\t/,$headers[$i]);
    my $gene_id = $newelement[0];
    $hash1{$gene_id} = [0,0,0,0,0,0,0,0,0,0,0,0,0]; # AMP-Binding, Acy_transf_1, Condensation, Ketoacyl-synt_C(BKS-C)
                                                    # PP-Binding, ketoacyl-synt(BKS-N), NAD_4, adh_short, DMAT
    $hash1{$gene_id}[13] = $newelement[1]; # gene_function
    $hash1{$gene_id}[14] = $newelement[2]; # contig/chromosome
    $hash1{$gene_id}[15] = $newelement[3]; # gene_order
    $hash1{$gene_id}[16] = $newelement[4]; # gene_start5
    $hash1{$gene_id}[17] = $newelement[5]; # gene_end3
    $hash1{$gene_id}[18] = $sequences[$i]; # protein sequences
}# End for loop

# Read in .htab files from $tempdir directory
my @files = <$tempdir/*.out.htab>;

# Loop through each HMMoutput.htab file and get the necessary information
foreach $filename1(@files){
    my @r = []; # store information from each file
    my $id = ''; # store gene_id
    my %hash2; # store domain content
    #print OUTFILE "$filename1\n"; # to determine the order of the domains i.e. order of columns in table.
    #print "$filename1\n";
    # Read in HMMouput file (from htab.pl), create hash2
    local $/="\n";  # going back to default setting '\n' being the seperator
    open(INFILE, $filename1) or die "Can't open $filename1";
    while (<INFILE>) {
	@r=split(/\t/,$_); # split the file by \t
	chomp @r; # remove extra \n
	$id = ">".$r[5]; # concatenate ">" infront of gene_id to match
        $hash2{$id}[0] = $r[14]; # total number of domains as read from htab.pl
	$hash2{$id}[1] = $r[11]; # domain score
    }#end while

    # Loop to check if gene_ids are the same in both files (inputFASTA file & HMMoutput)
    # & if the same add the info. from the first file to that
    # gene_id--->to prevent mismatches
    # create hash1 (hash containing all the information)
    foreach $k(keys %hash1){
	if(exists $hash2{$k}){
	    if($filename1 eq "$tempdir/hmm_Acyl_transf.out.htab"){
		$hash1{$k}[0] = $hash2{$k}[0];
	    }
	    if($filename1 eq "$tempdir/hmm_AMP.out.htab"){
                $hash1{$k}[1] = $hash2{$k}[0];
	    }
	    if($filename1 eq "$tempdir/hmm_Condensation.out.htab"){
                $hash1{$k}[2] = $hash2{$k}[0];
	    }
	    if($filename1 eq "$tempdir/hmm_Ketoacyl-synt_C.out.htab"){
                $hash1{$k}[3] = $hash2{$k}[0];
		$hash1{$k}[21] = $hash2{$k}[1];
	    }
	    if($filename1 eq "$tempdir/hmm_ketoacyl-synt.out.htab"){
                $hash1{$k}[4] = $hash2{$k}[0];
		$hash1{$k}[22] = $hash2{$k}[1];
	    }
	    if($filename1 eq "$tempdir/hmm_PP-binding.out.htab"){
                $hash1{$k}[5] = $hash2{$k}[0];
	    }
	    #if($filename1 eq "$tempdir/hmm_Thioesterase.out.htab"){
                #$hash1{$k}[6] = $hash2{$k}[0];
	    #}
	    if($filename1 eq "$tempdir/hmm_NAD_binding_4.out.htab"){
                $hash1{$k}[6] = $hash2{$k}[0];
	    }
	    if($filename1 eq "$tempdir/hmm_adh_short.out.htab"){
                $hash1{$k}[7] = $hash2{$k}[0];
	    }
	    if($filename1 eq "$tempdir/hmm_DMAT.out.htab"){
                $hash1{$k}[8] = $hash2{$k}[0];
	    }
	    if($filename1 eq "$tempdir/hmm_arom_pren_DMATS.out.htab"){
                $hash1{$k}[9] = $hash2{$k}[0];
	    }
	    if($filename1 eq "$tempdir/hmm_trp_dimet_allyl.out.htab"){
                $hash1{$k}[10] = $hash2{$k}[0];
            }
	    if($filename1 eq "$tempdir/hmm_Lys2.out.htab"){
                $hash1{$k}[11] = $hash2{$k}[0];
	    }
	    if($filename1 eq "$tempdir/hmm_alpha_am_amid.out.htab"){
                $hash1{$k}[12] = $hash2{$k}[0];
		$hash1{$k}[20] = $hash2{$k}[1];
	    }
	    #if($filename1 eq "$tempdir/hmm_Epimerase.out.htab"){
                #$hash1{$k}[9] = $hash2{$k}[0]
	    #}
	    #if($filename1 eq "$tempdir/hmm_KR.out.htab"){
                #$hash1{$k}[10] = $hash2{$k}[0];
	    #}

	    if($hash1{$k}[1]>0 && $hash1{$k}[2]>0 && $hash1{$k}[5]>0 && $hash1{$k}[0]>0 && $hash1{$k}[3]>0 && $hash1{$k}[4]>0){
                $hash1{$k}[19] = "HYBRID";
            }# end if loop
	    elsif($hash1{$k}[1]>0 && $hash1{$k}[2]>0 && $hash1{$k}[5]>0){ # determine if NRPS
		$hash1{$k}[19] = "NRPS";
	    }# end elsif loop
	    elsif($hash1{$k}[0]>0 && $hash1{$k}[3]>0 && $hash1{$k}[4]>0){ # determine if PKS
		$hash1{$k}[19] = "PKS";
	    }# end elsif
	    elsif(($hash1{$k}[1]>0 && $hash1{$k}[2]>0) || ($hash1{$k}[1]>0 && $hash1{$k}[5]>0) ||($hash1{$k}[2]>0 && $hash1{$k}[5]>0) || ($hash1{$k}[1]>0 && $hash1{$k}[6]>0) || ($hash1{$k}[1]>0 && $hash1{$k}[7]>0)){
		$hash1{$k}[19] = "NRPS-Like";
            }# end elsif
	    elsif(($hash1{$k}[0]>0 && $hash1{$k}[3]>0) || ($hash1{$k}[0]>0 && $hash1{$k}[4]>0) || ($hash1{$k}[3]>0 && $hash1{$k}[4]>0)){
		$hash1{$k}[19] = "PKS-Like";
	    }# end elsif
	    elsif($hash1{$k}[8]>0 && $hash1{$k}[9]>0 && $hash1{$k}[10]>0){
		$hash1{$k}[19] = "DMAT";
	    }# end elsif
	    else{ # neither NRPS nor PKS
		$hash1{$k}[19] = "NONE";
	    }# end else
	}# end if loop
    }# end foreach loop
    close (INFILE);
} # end foreach loop

print "Backbone_gene_id\tAnnotated_gene_function\tChromosome-Contig\tGene_order\t5'_end\t3'_end\tSMURF_backbone_gene_prediction\n";

# print the info. hash1
foreach my $keys(sort keys %hash1){
    #####FIX THIS CONDITION#########
    #####REMOVE EXTRA DOMAINS######ALSO FROM testscript2.pl
    if($hash1{$keys}[0]==0 && $hash1{$keys}[1]==0 && $hash1{$keys}[2]==0 && $hash1{$keys}[3]==0 && $hash1{$keys}[4]==0 && $hash1{$keys}[5]==0 && $hash1{$keys}[6]==0 && $hash1{$keys}[7]==0 && $hash1{$keys}[8]==0 && $hash1{$keys}[9]==0 && $hash1{$keys}[10]==0 && $hash1{$keys}[11]==0 && $hash1{$keys}[12]==0){ # loop to eliminate columns for which all domains are 0
   }# end if
    else{
	# eliminate "alpha" domain for NRPS-like enzyme using domain score cut-off value
	if(not defined $hash1{$keys}[20]){
            $hash1{$keys}[20] = 0;
        }# end if
        else{
            $hash1{$keys}[20] = $hash1{$keys}[20];
        }# end else
	if($hash1{$keys}[20]>0){
	    $hash1{$keys}[19] = "NONE";
	}# end if

	# eliminate "false positives" for PKS-like enzyme using domain score cut-off value
	if(not defined $hash1{$keys}[21]){
            $hash1{$keys}[21] = 0;
        }# end if
        else{
            $hash1{$keys}[21] = $hash1{$keys}[21];
        }# end else
	if(not defined $hash1{$keys}[22]){
            $hash1{$keys}[22] = 0;
        }# end if
        else{
            $hash1{$keys}[22] = $hash1{$keys}[22];
        }# end else
	if($hash1{$keys}[21]<0 || $hash1{$keys}[22]<0){
            $hash1{$keys}[19] = "NONE";
        }# end if

        # {gene_id},{AMP-binding},{Acyl_transf_1},{Condensation},{Ketoacyl-synt_C},{PP-binding},{Thioester-redct},
	# {Thioesterase},{ketoacyl-synt}
	#print OUTFILE "$keys\t$hash1{$keys}[0]\t$hash1{$keys}[1]\t$hash1{$keys}[2]\t$hash1{$keys}[3]\t$hash1{$keys}[4]\t$hash1{$keys}[5]\t$hash1{$keys}[6]\t$hash1{$keys}[7]\t$hash1{$keys}[8]\t$hash1{$keys}[9]\t$hash1{$keys}[10]\t$hash1{$keys}[11]\t$hash1{$keys}[12]\t$hash1{$keys}[13]\t$hash1{$keys}[19]\n";
	   if($hash1{$keys}[19] eq "NONE"){ # eliminate neither NRPS nor PKS
	   }#end if
	   else{ # print all NRPS's, PKS's
	       #print OUTFILE "$keys\t$hash1{$keys}[13]\t$hash1{$keys}[14]\t$hash1{$keys}[15]\t$hash1{$keys}[16]\t$hash1{$keys}[17]\t$hash1{$keys}[19]\n$hash1{$keys}[18]\n";
	       #print "$keys\t$hash1{$keys}[13]\t$hash1{$keys}[14]\t$hash1{$keys}[15]\t$hash1{$keys}[16]\t$hash1{$keys}[17]\t$hash1{$keys}[19]\n$hash1{$keys}[18]\n";
               #print OUTFILE "$keys\t$hash1{$keys}[13]\t$hash1{$keys}[14]\t$hash1{$keys}[15]\t$hash1{$keys}[16]\t$hash1{$keys}[17]\t$hash1{$keys}[19]\n";
	       $keys =~ m/\>(.*)/;
	       print "$1\t$hash1{$keys}[13]\t$hash1{$keys}[14]\t$hash1{$keys}[15]\t$hash1{$keys}[16]\t$hash1{$keys}[17]\t$hash1{$keys}[19]\n";
	   }#end else
        }# end else
}# end foreach
#close OUTFILE;

# Potential NRPS's have:
# 1. AMP-binding enzymes
# 2. Condenstation domain
# 3. Phosphopantetheine attachment site
# 4. Thioester reductase
# 5. Thioesterase domain
# 6. Acyl transferase domain
# 7. Starter unit: ACP transacylase (SAT) pmid: 17071746 (need to build hmm model)

# Potential PKS's have:
# 1. Phosphopantetheine attachment site
# 2. Beta-ketoacyl synthase, C-terminal domain (BKS-C)
# 3. Beta-ketoacyl synthase, N-terminal domain (BKS-N)

exit;











































