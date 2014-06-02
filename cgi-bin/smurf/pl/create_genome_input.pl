#!/usr/local/bin/perl -w

# Author: Fayaz Seifuddin
# Date Created: 05/17/2007
# Last Date Modified: 05/22/2007
# File:script1.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: Read in peptidefile & coordinatefile, create input FASTA file, STEP1 
# Usage: perl program.pl peptidefile coordinatefile

# Enforces rules about variables, references
# and subroutines. This requires all variables be declared as
# "my variables"
use strict;

# Varible declaration
my $filename;  	# Name of the input file1
my $filename1; 	# Name of input file2 
my $line;      	# Number of lines in input file1
my $line1;	# Number of lines in input file2
my @headers; 	# Array to store header list for file1
my @headers1; 	# Array to store header list for file2
my @sequences; 	# Array to store sequence list for file1
my @newelement; # Array to store new split headers for file2 
my @newelement1; # Array to store new split headers for file1
my $counter=-1; # Counter to go through the entire FASTA file1
my $counter1=-1; # Counter to go through the entire header/coordinate file2
my %hash1;	# Hash to hold gene_ids-keys & coordinate information-values from file1
my %hash2;	# Hash to hold gene_ids-keys & coordinate information-values from file2
my $k;		# Variable for gene_ids
my $j=0;        # Variable for chromosome/super-contig genes order
my $contig = -1; # Test variable to check if we're on the same contig
my $chromo;     # Variable to store chromosome number

# Read from a FASTA file containing multiple sequences
# Read in header file containing coordinates
# Open the files given as arguments on the command line
$filename = $ARGV[0];
$filename1 = $ARGV[1];

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

open(INFILE, $filename1) or die "Can't open $filename1";
# Read in coordinates from file2
while ($line1 = <INFILE>) {
  chomp $line1; 			    # Remove extra end line \n from EOF
  # If loop to skip to the next line if the first line in file is blank
  if ($line1 =~ /^\s*$/){
    next;
  }
  else{               
    $counter1++;                        # Increament counter to go to next line
    $headers1[$counter1]=($line1);      # Save header line
  }
}# End while loop
close(INFILE);

# Loop to create hash from inputfile1
for (my $i=0; $i<=scalar(@headers)-1; $i++){
  @newelement1 = split(/ /,$headers[$i]);
  my $gene_id1 = $newelement1[0];
  $hash1{$gene_id1} = $sequences[$i];	    
}# End for loop

# Loop to create hash from inputfile2
for (my $i=0; $i<=scalar(@headers1)-1; $i++){
  @newelement = split(/\t/,$headers1[$i]); # Split on | or Split on \t
  #print "$newelement[0]\t$newelement[1]\t$newelement[2]\t$newelement[3]\t$newelement[4]\n";
  my $gene_id = ">".$newelement[0]; # GeneID 
  my $gene_name = $newelement[4];  # Gene_name  
  my $chromosome_number = $newelement[1];  # Chromosome_number
  my $gene_start5 = $newelement[2];  # 5'end (Gene_start)
  my $gene_start3 = $newelement[3];  # 3'end (Gene_stop)
  # Add above values to hash as arrays
  $hash2{$gene_id}[0] = $gene_name;
  $hash2{$gene_id}[1] = $chromosome_number;
  $hash2{$gene_id}[2] = $gene_start5;
  $hash2{$gene_id}[3] = $gene_start3;	
}# End for loop

# Loop to check if gene_ids are the same in both files
# & if the same add the sequence from the first file to that
# gene_id--->to prevent mismatches
foreach $k(sort keys %hash2){
  if(exists $hash1{$k}){
    push(@{$hash2{$k}[4]},$hash1{$k});
  }# End if loop
}# End foreach loop

#open OUTFILE,"> OUT" or die "Can't open: OUT\n";

# Loop to print the new FASTA file/Sort by chromosome/contig, gene_start5,gene_id
foreach $k(sort {$hash2{$a}[1]<=>$hash2{$b}[1]}  sort {$hash2{$a}[2]<=>$hash2{$b}[2]} sort keys(%hash2)){
  $chromo = $hash2{$k}[1]; # Chromosome numbers
  if($contig==$chromo){  # Check if chromosome number is the same
    $j=$j+1;  # Increment Gene_Order if chromosome number is the same
    $hash2{$k}[5]=$j; # Add Gene_Order value to hash
  }#endif
  else{ # Re-set variable/gene_order to 1 if chromosome number changes and change contig number as well
    $j=1;
    $contig = $chromo;
    $hash2{$k}[5]=$j ;
  }#endelse
  if(not defined @{$hash2{$k}[4]}){
  }
  else{
    #print OUTFILE "$k\t$hash2{$k}[0]\t$hash2{$k}[1]\t$hash2{$k}[5]\t$hash2{$k}[2]\t$hash2{$k}[3]\n@{$hash2{$k}[4]}\n";
    print "$k\t$hash2{$k}[0]\t$hash2{$k}[1]\t$hash2{$k}[5]\t$hash2{$k}[2]\t$hash2{$k}[3]\n@{$hash2{$k}[4]}\n";
  }
}# End foreach loop

#close OUTFILE;
exit;
