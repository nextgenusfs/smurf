#!/usr/local/bin/perl -w

# Author: Fayaz Seifuddin, Nora Khaldi
# Date Created: 08/17/2007
# Last Date Modified: 08/17/2007
# File:sortclusterproblem.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: sort cluster problems for sorting according to gene_order
# Usage: perl sortclusterproblem.pl CLUSTER-DOMAIN-CONTENTOUT

use strict;

# Variable declarations
my @h;
my @s;
my @e_ids;
my @e_names;
my @clustenzyme_ids;
my @clustg_ids;
my @positions;
my @contigs;
my @gene_order;
my @starts_5;
my @ends_3;
my @clustg_dists;
my @clustg_functions;
my $domain;
my $domaincheck;
my @domains;
my @domainchecks;
my @newelements;
my %clusters;
my %clustersn;

# Open the files given as arguments on the command line
my $filename = $ARGV[0]; # CLUSTER-DOMAIN-CONTENT-OUT

&open_file($filename,\@h,\@s); # CLUSTER-DOMAIN-CONTENT-OUT file (call subroutine open_file)
# Loop to create ARRAYS from CLUSTER-DOMAIN-CONTENT-OUT file
for (my $j=0; $j<=scalar(@h)-1; $j++){
    $domain = '';
    $domaincheck = '';
    @newelements = split(/\t/,$h[$j]); # split each line by tab
    push(@clustenzyme_ids, $newelements[0]);
    push(@clustg_ids, $newelements[1]);
    push(@positions, $newelements[2]);
    push(@contigs, $newelements[3]);
    push(@gene_order, $newelements[4]);
    push(@starts_5, $newelements[5]);
    push(@ends_3, $newelements[6]);
    push(@clustg_dists, $newelements[7]);
    for (my $l=8; $l<scalar(@newelements); $l++){#changed 08.07.2008
        $domain = $domain.$newelements[$l];
        if($domain != 0){
            $domaincheck = 1;
        }
	else{
            $domaincheck = 0;
        }
    }
    push(@domains,$domain);
    push(@domainchecks, $domaincheck);
    push(@clustg_functions, $newelements[39]);
    $clusters{$clustenzyme_ids[$j]}{$gene_order[$j]}{$clustg_ids[$j]}{$positions[$j]}{$contigs[$j]}{$starts_5[$j]}{$ends_3[$j]}{$clustg_dists[$j]}{$domainchecks[$j]}{$clustg_functions[$j]} = 0;
}# End for loop

foreach my $c_e(sort keys %clusters){
    foreach my $gene_order(sort {$a <=> $b} keys %{$clusters{$c_e}}){
	foreach my $g_id(sort keys %{$clusters{$c_e}{$gene_order}}){
	    foreach my $pos(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}}){
		foreach my $chromos(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}}){
		    foreach my $start_5(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}}){
			foreach my $end_3(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}}){
			    foreach my $g_distance(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}}){
				foreach my $domain_score(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}{$g_distance}}){
				    foreach my $gene_function(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}{$g_distance}{$domain_score}}){
					print "$c_e\t$g_id\t$pos\t$chromos\t$gene_order\t$start_5\t$end_3\t$g_distance\t$domain_score\t$gene_function\n";
				    }#endforeeach10
				}#endforeach9
			    }#endforeach8
			}#endforeach7
		    }#endforeach6
		}#endforeach5
	    }#endforeach4
        }#endforeach3
    }#endforeach2
}#endforeach1

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
