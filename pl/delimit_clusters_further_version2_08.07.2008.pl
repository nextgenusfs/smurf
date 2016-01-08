#!/usr/bin/env perl -w

# Author: Fayaz Seifuddin, Nora Khaldi
# Date Created: 08/07/2008
# Last Date Modified: 08/07/2008
# File: new_testclusterstep4.2.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: Read delimited clusters, remove those genes with domainscore==0 from the tails/ends, find 
# overlapping clusters, merge overlapping clusters into one, assign cluster ID's, output new cluster file 
# Usage: perl new_testclusterstep4.2.pl DELIMITED CLUSTERS DOMAINCONTENTOUT-(backbone enzymes)

use strict;

# Variable declarations
my @h;
my @s;
my @newelements;
my @clustenzyme_ids;
my @clustg_ids;
my @positions;
my @contigs;
my @gene_order;
my @starts_5;
my @ends_3;
my @clustg_dists;
my @domainscores;
my @clustg_functions;
my %clusters;
my %clustersn;
my $domain;
my $bool = 'false';
my $domain_position;

###
my $offset = 0;
my $offset_flag = 'true';
my $flag = 'true';
###

# Open the files given as arguments on the command line
my $filename = $ARGV[0];  # DELIMITED CLUSTERS
my $filename1 = $ARGV[1];  # DOMAINCONTENTOUT

&open_file($filename,\@h,\@s); # CLUSTER-DOMAIN-CONTENT-OUT file (call subroutine open_file)
# Loop to create ARRAYS from CLUSTER-DOMAIN-CONTENT-OUT file
for (my $j=0; $j<=scalar(@h)-1; $j++){
    @newelements = split(/\t/,$h[$j]); # split each line by tab
    $newelements[0] =~ m/\>(.*)/;
    $newelements[0] = $1;
    push(@clustenzyme_ids, $newelements[0]);
    $newelements[1] =~ m/\>(.*)/;
    $newelements[1] = $1;
    push(@clustg_ids, $newelements[1]);
    push(@positions, $newelements[2]);
    push(@contigs, $newelements[3]);
    push(@gene_order, $newelements[4]);
    push(@starts_5, $newelements[5]);
    push(@ends_3, $newelements[6]);
    push(@clustg_dists, $newelements[7]);
    push(@domainscores, $newelements[8]);
    push(@clustg_functions, $newelements[9]);
    $clustersn{$clustenzyme_ids[$j]}{$gene_order[$j]}{$clustg_ids[$j]}{$positions[$j]}{$contigs[$j]}{$starts_5[$j]}{$ends_3[$j]}{$clustg_dists[$j]}{$domainscores[$j]}{$clustg_functions[$j]}=0;

}# End for loop

@h = [];
@s = [];

&open_file($filename1,\@h,\@s); # DOMAINCONTENTOUT file (call subroutine open_file)
# Loop to create ARRAYS from DOMAINCONTENTOUT file                                                                                                                               for(my $j=0; $j<=scalar(@h)-1; $j++) {
    my @element = split(/\t/,$h[$j]); # split headers by tab                                                                                                                         my $e_id = $element[0]; # NRPS,PKS,DMAT ids                                                                                                                                      my $e_name = $element[6]; # type of enzyme,enzyme name                 
    push(@e_ids,$e_id); # store NRPS,PKS,DMAT ids                                                                                                                                    push(@e_names, $e_name); # store type of enzyme, enzyme name...                                                                                                             
}#End for loop 

for(my $i=0; $i<=scalar(@clustg_ids)-1; $i++) {


}#end for






##############################################################################################################################################################################


for(my $i=0; $i<=scalar(@clustenzyme_ids)-1; $i++){
    if($positions[$i] == 0){
	$offset = 0;
	$offset_flag = 'true';
	$flag = 'true';
	##############DOWNSTREAM#########################
	for(my $j=20; $j>=0; --$j){
	    if(not defined $clustenzyme_ids[$i-$j] || $positions[$i-$j] || $contigs[$i-$j] || $gene_order[$i-$j] || $starts_5[$i-$j] || $ends_3[$i-$j] || $domainscores[$i-$j] || $clustg_ids[$i-$j] || $clustg_dists[$i-$j] || $clustg_functions[$i-$j]){
		next;
	    }#end if
	    else{
		if($clustenzyme_ids[$i] eq $clustenzyme_ids[$i-$j]){
		    if($domainscores[$i] == 0){
			$domainscores[$i] = 1;
		    }
		    if($offset_flag eq 'true'){
			$offset = $i-$j;
			$offset_flag = 'false';
			#print "OFFSET:", "\t", $offset, "\n";
		    }#end if
		    if($domainscores[$i-$j] == 1 && $flag eq 'true'){
			my $y = $i-$j;
			#print "FIRST ONE FOUND:", "\t", $y, "\n";
			if(($y - $offset) >= 4 && $positions[$i-$j] == 0){
			    #print "HOW MANY ZEROS:", "\t", $y - $offset, "\n";
			    for(my $r=$y-2; $r<=$i; $r++){
				$clusters{$clustenzyme_ids[$r]}{$gene_order[$r]}{$clustg_ids[$r]}{$positions[$r]}{$contigs[$r]}{$starts_5[$r]}{$ends_3[$r]}{$clustg_dists[$r]}{$domainscores[$r]}{$clustg_functions[$r]}=0;
				$clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
                            #print "$clustenzyme_ids[$r]\n";
                            }#end for
			}#end if
			    elsif(($y - $offset >= 10)){
                                #print "I AM IN THIS LOOP\n";
                                for(my $f=$y-2; $f<=$i; $f++){
                                    $clusters{$clustenzyme_ids[$f]}{$gene_order[$f]}{$clustg_ids[$f]}{$positions[$f]}{$contigs[$f]}{$starts_5[$f]}{$ends_3[$f]}{$clustg_dists[$f]}{$domainscores[$f]}{$clustg_functions[$f]}=0;
                                    $clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
                            #print "$clustenzyme_ids[$r]\n";
                                }#end for
			    }#end elsif
			elsif(($y - $offset) >= 4){
			    #print "HOW MANY ZEROS:", "\t", $y - $offset, "\n";
                            for(my $l=$y; $l<=$i; $l++){
                                $clusters{$clustenzyme_ids[$l]}{$gene_order[$l]}{$clustg_ids[$l]}{$positions[$l]}{$contigs[$l]}{$starts_5[$l]}{$ends_3[$l]}{$clustg_dists[$l]}{$domainscores[$l]}{$clustg_functions[$l]}=0;
                                $clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
                            #print "$clustenzyme_ids[$r]\n";
                            }#end for
			}#end elsif
			    #elsif(($y - $offset >= 10)){
				#print "I AM IN THIS LOOP\n";
				#for(my $f=$y-2; $f<=$i; $f++){
				    #$clusters{$clustenzyme_ids[$f]}{$gene_order[$f]}{$clustg_ids[$f]}{$positions[$f]}{$contigs[$f]}{$starts_5[$f]}{$ends_3[$f]}{$clustg_dists[$f]}{$domainscores[$f]}{$clustg_functions[$f]}=0;
				    #$clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
                            #print "$clustenzyme_ids[$r]\n";
				#}#end for
			    #}#end elsif
			else{
			    #print "LESS ZEROS:", "\t", $y - $offset, "\n";
			    for(my $t=$offset; $t<=$i; $t++){
			    $clusters{$clustenzyme_ids[$t]}{$gene_order[$t]}{$clustg_ids[$t]}{$positions[$t]}{$contigs[$t]}{$starts_5[$t]}{$ends_3[$t]}{$clustg_dists[$t]}{$domainscores[$t]}{$clustg_functions[$t]}=0;
			    $clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
			    #print "$clustenzyme_ids[$t]\n";
			    }#end for
			}#end else
			    $flag = 'false';
		    }#end if
		}#end if
	    }#end else
	}#end for
	##################UPSTREAM#####################
	$offset = 0;
        $offset_flag = 'true';
        $flag = 'true';
	for(my $j=20; $j>=0; --$j){
            if(not defined $clustenzyme_ids[$i+$j] || $positions[$i+$j] || $contigs[$i+$j] || $gene_order[$i+$j] || $starts_5[$i+$j] || $ends_3[$i+$j] || $domainscores[$i+$j] ||$clustg_ids[$i+$j] || $clustg_dists[$i+$j] || $clustg_functions[$i+$j]){
                next;
            }#end if
		else{
		    if($clustenzyme_ids[$i] eq $clustenzyme_ids[$i+$j]){
			if($domainscores[$i] == 0){
			    $domainscores[$i] = 1;
			}
			if($offset_flag eq 'true'){
			    $offset = $i+$j;
			    $offset_flag = 'false';
			    #print "OFFSET:", "\t", $offset, "\n";
			}#end if
			    if($domainscores[$i+$j] == 1 && $flag eq 'true'){
				my $y = $i+$j;
				#print "FIRST ONE FOUND:", "\t", $y, "\n";
				if(($offset - $y) >= 4 && $positions[$i+$j]==0){
				    #print "HOW MANY ZEROS:", "\t", $offset - $y, "\n";
				    for(my $r=$y+2; $r>=$i; $r--){
					$clusters{$clustenzyme_ids[$r]}{$gene_order[$r]}{$clustg_ids[$r]}{$positions[$r]}{$contigs[$r]}{$starts_5[$r]}{$ends_3[$r]}{$clustg_dists[$r]}{$domainscores[$r]}{$clustg_functions[$r]}=0;
					$clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
					 #print "$clustenzyme_ids[$r]\n";
				    }#end for
				}#end if
			     elsif(($offset - $y >= 10)){
				 #print "I AM HERE IN THIS LOOP\n";
                                for(my $g=$y+2; $g>=$i; $g--){
                                    $clusters{$clustenzyme_ids[$g]}{$gene_order[$g]}{$clustg_ids[$g]}{$positions[$g]}{$contigs[$g]}{$starts_5[$g]}{$ends_3[$g]}{$clustg_dists[$g]}{$domainscores[$g]}{$clustg_functions[$g]}=0;
                                    $clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
                            #print "$clustenzyme_ids[$r]\n";
                                }#end for
                            }#end elsif	    
			elsif(($offset - $y) >= 4){
                            #print "HOW MANY ZEROS:", "\t", $offset - $y, "\n";
                            for(my $l=$y; $l>=$i; $l--){
                                $clusters{$clustenzyme_ids[$l]}{$gene_order[$l]}{$clustg_ids[$l]}{$positions[$l]}{$contigs[$l]}{$starts_5[$l]}{$ends_3[$l]}{$clustg_dists[$l]}{$domainscores[$l]}{$clustg_functions[$l]}=0;
                                $clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
                            #print "$clustenzyme_ids[$r]\n";
                            }#end for
                        }#end elsif    
			    #elsif(($offset - $y >= 10)){
                                #for(my $g=$y-2; $g>=$i; $g--){
                                    #$clusters{$clustenzyme_ids[$g]}{$gene_order[$g]}{$clustg_ids[$g]}{$positions[$g]}{$contigs[$g]}{$starts_5[$g]}{$ends_3[$g]}{$clustg_dists[$g]}{$domainscores[$g]}{$clustg_functions[$g]}=0;
				    #$clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
                            #print "$clustenzyme_ids[$r]\n";
                                #}#end for
			    #}#end elsif
			        else{
				    #print "LESS ZEROS:", "\t", $offset - $y, "\n";
				    for(my $t=$offset; $t>=$i; $t--){
					$clusters{$clustenzyme_ids[$t]}{$gene_order[$t]}{$clustg_ids[$t]}{$positions[$t]}{$contigs[$t]}{$starts_5[$t]}{$ends_3[$t]}{$clustg_dists[$t]}{$domainscores[$t]}{$clustg_functions[$t]}=0;
				       $clusters{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainscores[$i]}{$clustg_functions[$i]}=0;
					#print "$clustenzyme_ids[$t]\n";
				    }#end for
				}#end else
				$flag = 'false';
			    }#end if
		      }#end if
		 }#end else
	  }#end for
    }#end if
}#end for

#foreach my $c_e(sort keys %clusters){
    #foreach my $c_g(sort keys %{$clusters{$c_e}}){
	#foreach my $pos(sort keys %{$clusters{$c_e}{$c_g}}){
	    #foreach my $dis(sort keys %{$clusters{$c_e}{$c_g}{$pos}}){
		#foreach my $dom(sort keys %{$clusters{$c_e}{$c_g}{$pos}{$dis}}){
		    #foreach my $bool(sort keys %{$clusters{$c_e}{$c_g}{$pos}{$dis}{$dom}}){
			#print "$c_e\t$c_g\t$pos\t$dis\t$dom\n";
		    #}#end foreach
		#}#end foreach5
	    #}#end foreach4
	#}#end foreach3
    #}#end foreach2
#}#end foreach1

foreach my $c_e(sort keys %clusters){
    foreach my $c_g(sort keys %{$clusters{$c_e}}){
	foreach my $pos(sort keys %{$clusters{$c_e}{$c_g}}){
	    foreach my $contigs(sort keys %{$clusters{$c_e}{$c_g}{$pos}}){
		foreach my $gene_order(sort keys %{$clusters{$c_e}{$c_g}{$pos}{$contigs}}){
		    foreach my $start_5(sort keys %{$clusters{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}}){
			foreach my $end_3(sort keys %{$clusters{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}{$start_5}}){
			    foreach my $dis(sort keys %{$clusters{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}{$start_5}{$end_3}}){
				foreach my $dom(sort keys %{$clusters{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}{$start_5}{$end_3}{$dis}}){
				    foreach my $g_function(sort keys %{$clusters{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}{$start_5}{$end_3}{$dis}{$dom}}){
					    #print "$c_e\t$c_g\t$pos\t$contigs\t$gene_order\t$start_5\t$end_3\t$dis\t$dom\t$g_function\n";
					}#end foreach10
				}#end foreach9
			    }#end foreach8
			}#end foreach7
		    }#end foreach6
		}#end foreach5
	    }#end foreach4
	}#end foreach3
    }#end foreach2
}#end foreach1

    my $i=1;
    foreach my $c_e(sort keys %clusters){
	print "Cluster:" ,$i, "\n";
	print "Backbone_gene_id\tGene_id\tGene_positions\tChromosome-Contig\tGene_order\t5'end\t3'end\tGene_distance\tDomain_score\tAnnotated_gene_function\n";
	$i++;
	#$c_e =~ m/\>(.*)/;
	#$c_e = $1;
	foreach my $gene_order(sort {$a <=> $b} keys %{$clusters{$c_e}}){
	    foreach my $g_id(sort keys %{$clusters{$c_e}{$gene_order}}){
		foreach my $pos(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}}){
		    foreach my $chromos(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}}){
			foreach my $start_5(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}}){
			    foreach my $end_3(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}}){
				foreach my $g_distance(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}}){
				    foreach my $domain_score(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}{$g_distance}}){
					foreach my $gene_function(sort keys %{$clusters{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}{$g_distance}{$domain_score}}){
					    #$c_e =~ m/\>(.*)/;
					    #$c_e = $1;

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
