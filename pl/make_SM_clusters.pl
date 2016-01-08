#!/usr/bin/env perl -w

# Author: Fayaz Seifuddin, Nora Khaldi
# Date Created: 06/11/2007
# Last Date Modified: 06/11/07
# File:script5.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: Read in inputFASTA file and DOMAINCONTENTOUT (NRPS,PKS,DMAT) files to get 20 genes upstream and 20 genes downstream from the core, calculate GENE DISTANCES from each other i.e. NRPS,PKS,DMAT, STEP3
# +/- 20 genes from NRPS,PKS,DMAT
# Usage: perl script5.pl inputFASTA DOMAINCONTENTOUT

use strict;

# Variable declarations
my %hash; # store info. from inputFASTA file
my %hash1; # store info. from DOMAINCONTENT file
my %hash2; # store info. for +/- 20 genes from NRPS,PKS,DMAT
my %hash3; # store info. for everything including position of NRPS,PKS,DMAT
my @h; # store info. for headers from inputFASTA file
my @s; # store info. for sequences from inputFASTA file
my @ids; # store info. for gene_ids from %hash (inputFASTA file)-sorted by chromosome/contig,gene_order,gene_start5,gene_id
my @ids1; # store info. for gene_ids from %hash1 (DOMAINCONTENTOUT)-NRPS,PKS,DMAT ids
my @chromosomes; # store info. for chromosomes from %hash (inputFASTA file)-sorted by chromosome/contig,gene_order,gene_start5,gene_id
my @gene_starts5; # store info. for gene_starts (5'starts) from %hash (inputFASTA file)-sorted by chromosome/contig,gene_order,gene_start5,gene_id
my @gene_ends3; # store info. for gene_ends (3'ends) from %hash (inputFASTA file)-sorted by chromosome/contig,gene_order,gene_start5,gene_id
my $gene_distance; # store info. for 5'-3'
my $gene_distance1; # store info. for 3'-5'
my $gene_distance2; # store info. for 3'-3'
my $gene_distance3; # store info. for 5'-5'
my @gene_distances; # store info. for all gene_distances above to calculate minimum gene_distance amongst all pairwise combinations

# Read from an inputFASTA file containing multiple sequences and headers
# Read from a DOMAINCONTENTOUT file containing predicted NRPS,PKS,DMAT
# Open the files given as arguments on the command line
my $filename = $ARGV[0];  # inputFASTA file
my $filename1 = $ARGV[1]; # DOMAINCONTENTOUT file

&open_file($filename,\@h,\@s); # inputFASTA file (call subroutine open_file)
# Loop to create hash from inputFASTA file
for (my $i=0; $i<=scalar(@h)-1; $i++){
    my @element = split(/\t/,$h[$i]); # split headers by tab
    my $gene_id = $element[0]; # gene_id
    my $gene_name = $element[1];  # gene_name/function
    my $chromosome_number = $element[2];  # chromosome_number/contig
    my $gene_order = $element[3]; # gene_order
    my $gene_start5 = $element[4];  # 5'end (gene_start)
    my $gene_start3 = $element[5];  # 3'end (gene_stop)
    # Add above values to hash as arrays
    $hash{$gene_id}[0] = $gene_name;
    $hash{$gene_id}[1] = $chromosome_number;
    $hash{$gene_id}[2] = $gene_order;
    $hash{$gene_id}[3] = $gene_start5;
    $hash{$gene_id}[4] = $gene_start3;
    $hash{$gene_id}[5] = $s[$i]; # sequence
}# End for loop

# Sort by chromosome/contig,gene_order,gene_start5,gene_id
# store gene_ids in an array in this order (order is GOOD here)
foreach my $k(sort {$hash{$a}[1]<=>$hash{$b}[1]} sort {$hash{$a}[2]<=>$hash{$b}[2]} sort {$hash{$a}[3]<=>$hash{$b}[3]} sort keys(%hash)){
    #print "$k\t$hash{$k}[1]\t$hash{$k}[2]\t$hash{$k}[3]\t$hash{$k}[4]\n";
    chomp $k; # remove extra \n
    push(@ids,$k); # store gene_ids
    push(@chromosomes, $hash{$k}[1]); # store chromosome_number
    push(@gene_starts5, $hash{$k}[3]); # store gene_starts (5'starts)
    push(@gene_ends3, $hash{$k}[4]); # store gene_ends (3'ends)
}# End foreach loop

# reset headers and sequences arrays to store info. from DOMAINCONTENTOUT file
@h = [];
@s = [];

&open_file($filename1,\@h,\@s); # DOMAINCONTENTOUT file (call subroutine open_file)
# Loop to create hash1 from DOMAINCONTENTOUT file
for (my $j=0; $j<=scalar(@h)-1; $j++){
    my @element1 = split(/\t/,$h[$j]); # split headers by tab
    my $gene_id1 = $element1[0]; # NRPS,PKS,DMAT ids
    #$hash1{$gene_id1} = $s[$j]; # assign sequence
    chomp $gene_id1; # remove extra \n
    push(@ids1,$gene_id1); # store NRPS,PKS,DMAT ids
}# End for loop

# compare @ids1(NRPS,PKS,DMAT ids) and @ids(gene_ids) to select +/- 20 genes from core
foreach my $a(@ids1){ # iterate through NRPS,PKS,DMAT ids
    for(my $b=0; $b<=scalar(@ids)-1; $b++){ # iterate through all gene_ids
	chomp $a; # remove extra \n
	if($a eq $ids[$b]){ # if NRPS,PKS,DMAT ids equal gene_ids
	    chomp $ids[$b]; # remove extra \n
	    for(my $c=1; $c<=20; $c++){ # select 20 genes downstream
		if(not defined $ids[$b-$c]){ # if at the beginning of the genome and not enough (<20) genes downstream add this gene_id
		    $c=21; # exit loop to select 20 genes downstream
		}#end if
		elsif($chromosomes[$b] != $chromosomes[$b-$c]){ # check if flanking genes are on the same chromosome for a cluster
		    $c=21;  # exit loop to select 20 genes downstream
		}#end if
		else{ # add flanking gene information to a hash
		    $gene_distance  = $gene_starts5[$b-$c] - $gene_ends3[$b-$c+1]; # 5'-3'
		    $gene_distance1 = $gene_ends3[$b-$c] - $gene_starts5[$b-$c+1]; # 3'-5'
		    $gene_distance2 = $gene_ends3[$b-$c] - $gene_ends3[$b-$c+1]; # 3'-3'
		    $gene_distance3 = $gene_starts5[$b-$c] - $gene_starts5[$b-$c+1]; #5'-5'
		    # append all gene_distances to an array to select minimum
		    @gene_distances = (abs($gene_distance),abs($gene_distance1),abs($gene_distance2),abs($gene_distance3));
		    # Call subroutine min_and_max to select minimum and maximum
		    my @ret = min_and_max(@gene_distances);
		    #print abs($gene_distance),"\t",abs($gene_distance1),"\t",abs($gene_distance2),"\t", abs($gene_distance3), "\t",$ret[0],"\n";
                    #print "$gene_starts5[$b-$c]\t$gene_ends3[$b-$c+1]\n";
		    $hash2{$a}{$ids[$b-$c]}{$c}{abs($ret[0])}=0; # add {NRPS,PKS,DMAT ids}{gene_ids}{gene_position from NRPS,PKS,DMAT}to hash2
		    #print $a,"\t",$ids[$b-$c],"\t",$c,"\t", abs($ret[0]), "\n";
	        }#end else
	    }#end for
	}#end if
    }#end for
}#end foreach

# compare @ids1(NRPS,PKS,DMAT ids) and @ids(gene_ids) to select +/- 20 genes from core
foreach my $a(@ids1){ # iterate through NRPS,PKS,DMAT ids
    for(my $b=0; $b<=scalar(@ids)-1; $b++){ # iterate through all gene_ids
	chomp $a; # remove extra \n
	if($a eq $ids[$b]){ # if NRPS,PKS,DMAT ids equal gene_ids
	    for(my $c=1; $c<=20; $c++){ # select 20 genes upstream
                if(not defined $ids[$b+$c]){ # if at the end of the genome and not enough (<20) genes upstream add this gene_id
		    $c=21; # exit loop to select 20 genes upstream
		}#end if
                elsif($chromosomes[$b] != $chromosomes[$b+$c]){ # check if flanking genes are on the same chromosome for a cluster
		    $c=21; # exit loop to select 20 genes upstream
		}#end if
                else{ # add flanking gene information to a hash
		    $gene_distance = $gene_starts5[$b+$c] - $gene_ends3[$b+$c-1]; # 5'-3'
		    $gene_distance1 = $gene_ends3[$b+$c] - $gene_starts5[$b+$c-1]; # 3'-5'
                    $gene_distance2 = $gene_ends3[$b+$c] - $gene_ends3[$b+$c-1]; # 3'-3'
                    $gene_distance3 = $gene_starts5[$b+$c] - $gene_starts5[$b+$c-1]; #5'-5'
		    # append all gene_distances to an array to select minimum
		    @gene_distances = (abs($gene_distance),abs($gene_distance1),abs($gene_distance2),abs($gene_distance3));
		    # Call subroutine min_and_max to select minimum and maximum
                    my @ret = min_and_max(@gene_distances);
                    #print $gene_distance,"\t",$gene_distance1,"\t",$gene_distance2,"\t",$gene_distance3, "\t",$ret[0],"\n";
		    #print "$gene_distance\n";
                    #print "$gene_starts5[$b+$c]\t$gene_ends3[$b+$c-1]\n";
		    $hash2{$a}{$ids[$b+$c]}{-$c}{abs($ret[0])}=0; # add {NRPS,PKS,DMAT ids}{gene_ids}{gene_position from NRPS,PKS,DMAT} to hash2
		    #print $a,"\t",$ids[$b+$c],"\t",-$c,"\t", abs($ret[0]), "\n";
		}#end else
            }#end for
	}#end if
    }#end for
}#end foreach

# add info. abt NRPS,PKS,DMAT at location '0' to a new hash (%hash3)
foreach my $enzyme(sort keys %hash2){
	foreach my $gid(sort keys %{$hash2{$enzyme}}){
		foreach my $loc(sort keys %{$hash2{$enzyme}{$gid}}){
		    foreach my $gdist(sort keys %{$hash2{$enzyme}{$gid}{$loc}}){
			#print "$enzyme\t$gid\t$loc\t$gdist\n";
			$hash3{$enzyme}{$enzyme}{0}{0} = $hash{$enzyme};
			$hash3{$enzyme}{$gid}{$loc}{$gdist} = $hash{$gid};
		    }#end foreach {$hash2{$enzyme}{$gid}{$loc}}
		}#end foreach {$hash2{$enzyme}{$gid}}
	}#end foreach {$hash2{$enzyme}}
}#end foreach {%hash2}

# write cluster-info to output file
#open OUTFILE,"> GENEDISTOUT" or die "Can't open: OUT\n";

# create hash3 and print cluster information to OUTFILE, GENEDISTOUT
foreach my $enzyme(sort keys %hash3){
	foreach my $gid(sort keys %{$hash3{$enzyme}}){
		foreach my $loc(sort keys %{$hash3{$enzyme}{$gid}}){
		    foreach my $gdist(sort keys %{$hash3{$enzyme}{$gid}{$loc}}){
			#print "$enzyme\t$gid\t$loc\n";
			#print OUTFILE "$enzyme\t$gid\t$loc\t$gdist\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[0]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[1]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[2]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[3]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[4]\n$hash3{$enzyme}{$gid}{$loc}{$gdist}[5]\n";
			print "$enzyme\t$gid\t$loc\t$gdist\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[0]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[1]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[2]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[3]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[4]\n$hash3{$enzyme}{$gid}{$loc}{$gdist}[5]\n";
                        #print OUTFILE "$enzyme\t$gid\t$loc\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[0]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[1]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[2]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[3]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[4]\n";
                        #print "$enzyme\t$gid\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[0]\t$hash3{$enzyme}{$gid}{$loc}{$gdist}[1]\n";
			#print "$enzyme\t$gid\t$loc\t$gdist\n";
		    }#end foreach {$hash2{$enzyme}{$gid}{$loc}}
		}#end foreach {$hash2{$enzyme}{$gid}}
	}#end foreach {$hash2{$enzyme}}
}#end foreach {%hash2}

#close OUTFILE;

exit;

########################################################################
#
#subroutine to open files, filenames provided as command line arguments
#
########################################################################
sub open_file{
    my($filename, $header_list, $sequence_list) = @_;
    my $line; # read one-line at time from <INFILE>
    my $counter = -1; # iterate through all headers and sequences in <INFILE>
    # open file
    open(INFILE, $filename) or die "Can't open $filename";
    while ($line = <INFILE>) {
	chomp $line;                          # Remove extra end line \n from EOF
    # If loop to skip to the next line if the first line in file is blank
	if ($line =~ /^\s*$/){
	    next;
	}#end if
	elsif ($line =~ /^>/) {                # Line starts with a ">"
	    $counter++;                        # Increament counter to go to next line
	    $$header_list[$counter]=($line);   # Save header line
	    $$sequence_list[$counter] = '';    # Start with an empty sequence
	}#end elsif
	else {
	    $$sequence_list[$counter]=$$sequence_list[$counter] . $line;     # Add line to end of sequence
	    $$sequence_list[$counter] =~ s/(\s|\n)//g;   # Remove all white space and end line (\n) from sequence
	}#end else
    }#end while
}#end sub open_file

########################################################################
#
#subroutine to select minimum and maximum gene_distances
#
########################################################################
sub min_and_max
{
    my (@numbers);
    @numbers = @_;
    my ($min, $max);
    $min = $numbers[0];
    $max = $numbers[0];
    foreach my $i (@numbers)
    {
        if($i > $max)
        {
            $max = $i;
        }
        elsif($i < $min)
        {
            $min = $i;
        }
    }
    return ($min, $max);
}#end subroutine min_and_max

