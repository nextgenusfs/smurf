#!/usr/local/bin/perl -w

# Author: Fayaz Seifuddin, Nora Khaldi
# Date Created: 08/17/2007
# Last Date Modified: 08/17/2007
# File:testclustersalldelimit.version2.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: Define delimited clusters for all clusters in a genome
# Usage: perl testclustersalldelimit.version2.pl DOMAINCONTENTOUT NEW SORTED BY GENE_ORDER CLUSTER-DOMAIN-CONTENTOUT

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
my @domainchecks;
my @newelements;
my %clusters;
my %clustersn;

my @distances;
my @max_min;
my $max_distance;

# Open the files given as arguments on the command line
my $filename = $ARGV[0];  # DOMAINCONTENTOUT
my $filename1 = $ARGV[1]; # CLUSTER-DOMAIN-CONTENT-OUT

&open_file($filename,\@h,\@s); # DOMAINCONTENTOUT file (call subroutine open_file)

# Loop to create ARRAYS from DOMAINCONTENTOUT file
for(my $j=0; $j<=scalar(@h)-1; $j++) {
  my @element = split(/\t/,$h[$j]); # split headers by tab
  my $e_id = $element[0]; # NRPS,PKS,DMAT ids
  my $e_name = $element[6]; # type of enzyme, enzyme name
  push(@e_ids,$e_id); # store NRPS,PKS,DMAT ids
  push(@e_names, $e_name); # store type of enzyme, enzyme name...
} # End for loop

@h = [];
@s = [];

&open_file($filename1,\@h,\@s); # CLUSTER-DOMAIN-CONTENT-OUT file (call subroutine open_file)

# Loop to create ARRAYS from CLUSTER-DOMAIN-CONTENT-OUT file
for(my $j=0; $j<=scalar(@h)-1; $j++) {
  @newelements = split(/\t/,$h[$j]); # split each line by tab
  push(@clustenzyme_ids, $newelements[0]);
  push(@clustg_ids, $newelements[1]);
  push(@positions, $newelements[2]);
  push(@contigs, $newelements[3]);
  push(@gene_order, $newelements[4]);
  push(@starts_5, $newelements[5]);
  push(@ends_3, $newelements[6]);
  push(@clustg_dists, $newelements[7]);
  push(@domainchecks, $newelements[8]);
  push(@clustg_functions, $newelements[9]);
  $clusters{$clustenzyme_ids[$j]}{$clustg_ids[$j]}{$positions[$j]}{$contigs[$j]}{$gene_order[$j]}{$starts_5[$j]}{$ends_3[$j]}{$clustg_dists[$j]}{$domainchecks[$j]}{$clustg_functions[$j]} = 0;
}# End for loop

for(my $i=0; $i<=scalar(@clustg_ids)-1; $i++) {
  for(my $k=0; $k<=scalar(@e_ids)-1; $k++) {
    ###ONLY NRPS NEXT TO ANY OTHER BACKBONE GENE IN THE SAME CLUSTER"###############
    if(($e_ids[$k] eq $clustg_ids[$i]) && ($e_names[$k] eq "NRPS" || $e_names[$k] eq "NRPS-Like") && ($positions[$i]>0 && $positions[$i]<=6) || (($e_ids[$k] eq $clustg_ids[$i]) && ($e_names[$k] eq "NRPS" || $e_names[$k] eq "NRPS-Like") && ($positions[$i]<0 && $positions[$i]>=-6))){
      my $var;
      my $cpt=0;

      $var = $positions[$i];

      if($var>=0){ 
        for(my $h=0; $h<=$var; $h++) {
          if(not defined $clustenzyme_ids[$i+$h]) {
            next;
          } # end if(not defined $clustenzyme_ids[$i+$h])
			$clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$h]}{$clustg_ids[$i+$h]}{$positions[$i+$h]}{$contigs[$i+$h]}{$starts_5[$i+$h]}{$ends_3[$i+$h]}{$clustg_dists[$i+$h]}{$domainchecks[$i+$h]}{$clustg_functions[$i+$h]}=0;
        } #end for(my $h=0; $h<=$var; $h++)

        for(my $t=1; $t<=20; $t++){
			 #print "$t\n";
			if(not defined $clustenzyme_ids[$i-$t]){
			    next;
			}
			 if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i-$t]){#IGNORE
			     #$t=21;
			 }#end if
			 elsif($positions[$i-$t]>=0){
			     #print "$positions[$i-$t]\n";
			     if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i-$t]){
				 $t=21;
				 #last;
			     }#end if
			     if($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]<=3814){
				 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
				 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
				 $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i-$t]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
			     	 #if($max_distance<$clustg_dists[$i-$t]){
				     #$max_distance=$clustg_dists[$i-$t];
				 #}#end if
			     }#end if
			     elsif($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]>=3814){
				 $t=21;
			     }#end elsif	 
		             elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]<=3814){
				 $cpt=$cpt+1;
				 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
				 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
				 $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i-$t]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
				 #print "$cpt\n";
				 #if($max_distance<$clustg_dists[$i-$t]){
				     #$max_distance=$clustg_dists[$i-$t];
				 #}#end if
			     }#end elsif
		             elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]>=3814){
				 $t=21;
			     }#end elsif
			     else{
			     }#end else;
		             if($cpt>=10){
				 $cpt=0; 
				 $t=21;
			     }#end if		  
			 }#end elsif
			 #elsif($positions[$i-$t]<=0){#this case should NOT happen because $position[$i-$t] ALWAYS POSITIVE
			     #print "$positions[$i+$t]\n";
			     #if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){
				 #$t=21;
			     #}#end if
                             #if($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]<2000){
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$gene_order[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
				 #if($max_distance<$clustg_dists[$i+$t]){
				     #$max_distance=$clustg_dists[$i+$t];
				 #}#end if
                             #}#end if
                             #elsif($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]>2000){
				 #$t=21;
			     #}#end elsif
                             #elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]<2000){
                                 #$cpt=$cpt+1;
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$gene_order[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
				 #print "$cpt\n";
                                 
				 #if($max_distance<$clustg_dists[$i+$t]){
				     #$max_distance=$clustg_dists[$i+$t];
				 #}#end if
                             #}#end elsif
                             #elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]>2000){
				 #$t=21;
			     #}#end elsif
                             #else{
			     #}#end else;
                             #if($cpt>=6){
				 #$cpt=0; 
				 #$t=21;
			     #}#end if
                         #}#end elsif   
		    }#end for
		    $cpt=0;
		    for(my $t=$var+1; $t<=20+$var; $t++){
                         #print "$t\n";
			if(not defined $clustenzyme_ids[$i+$t]){
			    next;
			}
			 if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){#IGNORE
			     #$t=21+$var;
			 }#end if
			 elsif($positions[$i+$t]<=0){
			     #print "$positions[$i+$t]\n";
			     if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){
				 $t=21+$var;
				 #last;
			     }#end if
			     if($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]<=3814){
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
				 $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$t]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
				 #if($max_distance<$clustg_dists[$i+$t]){
				     #$max_distance=$clustg_dists[$i+$t];
				 #}#end if
			     }#end if
			     elsif($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]>=3814){
				 $t=21+$var;
			     }#end elsif
			     elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]<=3814){
				 $cpt=$cpt+1;
				 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
				 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
				 $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$t]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
				 #print "$cpt\n";
				 #if($max_distance<$clustg_dists[$i+$t]){
				     #$max_distance=$clustg_dists[$i+$t];
				 #}#end if
			     }#end elsif
			     elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]>=3814){
				 $t=21+$var;
			     }#end elsif
			     else{
			     }#end else;	 
			     if($cpt>=10){
				 $cpt=0; 
				 $t=21+$var;
			     }#end id
                         }#end elsif
			 #elsif($positions[$i+$t]>=0){#this case should NOT happen because $position[$i+$t] ALWAYS NEGATIVE
			     #print "$positions[$i-$t]\n";
			     #if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i-$t]){
				 #$t=21+$var;
			     #}#end if
                             #if($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]<2000){
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$gene_order[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
				 #if($max_distance<$clustg_dists[$i-$t]){
				     #$max_distance=$clustg_dists[$i-$t];
				 #}#end if
                             #}#end if
                             #elsif($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]>2000){
				 #$t=21+$var;
			     #}#end elsif
                             #elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]<2000){
                                 #$cpt=$cpt+1;
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$gene_order[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
				 #print "$cpt\n";
                                 #if($max_distance<$clustg_dists[$i-$t]){
				     #$max_distance=$clustg_dists[$i-$t];
				 #}#end if
                             #}#end elsif
                             #elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]>2000){
				 #$t=21+$var;
			     #}#end elsif
			     #else{
			     #}#end else;
                             #if($cpt>=6){
				 #$cpt=0; 
				 #$t=21+$var;
			     #}#end if
                         #}#end elsif    
		    }#end for
		    $cpt=0;
	    }#end if
	    elsif($var<=0){
		$cpt = 0;
		#print "VAR:", $var, "\t", $clustenzyme_ids[$i], "\n";
		#print "$i\n";
		for(my $p=0; $p>=$var; --$p){
                    #print "$p\n";
		    #if($positions[$i+$p]>0){#SHOULD NOT HIT THIS CONDITION BECAUSE $positions[$i+$p] ALWAYS NEGATIVE
			#print "$positions[$i-$p]\n";
			#$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$p]}{$positions[$i-$p]}{$clustg_dists[$i-$p]}{$domainchecks[$i-$p]} = 0;
			#$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$p]}{$positions[$i-$p]}{$contigs[$i-$p]}{$gene_order[$i-$p]}{$starts_5[$i-$p]}{$ends_3[$i-$p]}{$clustg_dists[$i-$p]}{$domainchecks[$i-$p]}{$clustg_functions[$i-$p]}=0;
		    #}#end if
		    #else{
			#print "$positions[$i+$p]\n";
			#$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$p]}{$positions[$i+$p]}{$clustg_dists[$i+$p]}{$domainchecks[$i+$p]} = 0;
		    if(not defined $clustenzyme_ids[$i+$p]){
			next;
		    }
			$clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$p]}{$clustg_ids[$i+$p]}{$positions[$i+$p]}{$contigs[$i+$p]}{$starts_5[$i+$p]}{$ends_3[$i+$p]}{$clustg_dists[$i+$p]}{$domainchecks[$i+$p]}{$clustg_functions[$i+$p]}=0;
			#push(@distances,$clustg_dists[$i+$p]);
		    #}#end else
                }#end for
		       #@max_min = min_and_max(@distances);
		       #$max_distance = $max_min[1];
		       #print "$max_distance\n";
		for(my $t=1; $t<=20; $t++){
                        #print "$t\n";
		    if(not defined $clustenzyme_ids[$i+$t]){
			next;
		    }
		        if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){#IGNORE
			    #$t=21;
			}#end if
			elsif($positions[$i+$t]<=0){
			    #print "$positions[$i+$t]\n";
			    if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){
				$t=21;
			    }#end if
                            if($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]<=3814){
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                                $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$t]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0; 
				#if($max_distance<$clustg_dists[$i+$t]){
				     #$max_distance=$clustg_dists[$i+$t];
				#}#end if
                             }#end if
                             elsif($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]>=3814){
				 $t=21;
			     }#end elsif
                             elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]<=3814){
                                 $cpt=$cpt+1;
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                                 $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$t]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
                                 #print "$cpt\n";
                                 #if($max_distance<$clustg_dists[$i+$t]){
				     #$max_distance=$clustg_dists[$i+$t];
				 #}#end if
                             }#end elsif
                             elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]>=3814){
				 $t=21;
			     }#end elsif
                       if($cpt>=10){
				 $cpt=0; 
				 $t=21;
			     }#end if 
			}#end elsif
		      #elsif($positions[$i+$t]>=0){#this CONDITION WILL NEVER BE TRUE BECAUSE $positions[$i+$t] is ALWAYS NEGATIVE
			  #print "$positions[$i-$t]\n";
			    #if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i-$t]){
				#$t=21;
			    #}#end if
                            #if($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]<2000){
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                                  #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$gene_order[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
                                 #if($max_distance<$clustg_dists[$i-$t]){
				     #$max_distance=$clustg_dists[$i-$t];
				 #}#end if
                             #}#end if
                             #elsif($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]>2000){
				 #$t=21;
			     #}#end elsif
                             #elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]<2000){
                                 #$cpt=$cpt+1;
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                                  #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$gene_order[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
                                 #print "$cpt\n";
                                 #if($max_distance<$clustg_dists[$i-$t]){
				     #$max_distance=$clustg_dists[$i-$t];
				 #}#end if
                             #}#end elsif
                             #elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]>2000){
				 #$t=21;
			     #}#end elsif
                             #if($cpt>=6){
				 #$cpt=0; 
				 #$t=21;
			     #}#end if
                        #}#end elsif	    
	        }#end for
		$cpt=0;
		for(my $t=$var-1; $t>-21+$var; --$t){
                        #print "$t\n";
		    if(not defined $clustenzyme_ids[$i+$t]){
			next;
		    }
                        if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){#IGNORE
                            #$t=-22+$var;
                        }#end if
			elsif($positions[$i+$t]>=0){
			    #print "$positions[$i+$t]\n";
			    if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){
				$t=-22+$var;
			    }#end if
			    if($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]<=3814){
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                                  $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$t]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
                                 #if($max_distance<$clustg_dists[$i+$t]){
				     #$max_distance=$clustg_dists[$i+$t];
				 #}#end if
			    }#end if
			    elsif($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]>=3814){
				$t=-22+$var;
			    }#end elsif
			    elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]<=3814){
				 $cpt=$cpt+1;
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                                 $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$t]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
                                 #print "$cpt\n";
			         #if($max_distance<$clustg_dists[$i+$t]){
				     #$max_distance=$clustg_dists[$i+$t];
				 #}#end if
			    }#end elsif
			    elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]>=3814){
				$t=-22+$var;
			    }#end elsif
			    if($cpt>=10){
				$cpt=0; 
				$t=-22+$var;
			    }#end if
			}#end elsif
			#elsif($positions[$i+$t]<=0){#IT WILL NEVER HIT THIS CONDITION BECAUSE $positions[$i+$t] ALWAYS POSITIVE
			    #print "$positions[$i-$t]\n";
			    #if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i-$t]){
				#$t=-22+$var;
			    #}#end if
                            #if($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]<2000){
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                                  #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$gene_order[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
                                 #if($max_distance<$clustg_dists[$i-$t]){
				     #$max_distance=$clustg_dists[$i-$t];
				 #}#end if
                            #}#end if
                            #elsif($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]>2000){
				#$t=-22+$var;
			    #}#end elsif
                            #elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]<2000){
                                 #$cpt=$cpt+1;
                                 #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
                                 #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                                  #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$gene_order[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
                                 #print "$cpt\n";
                                 #if($max_distance<$clustg_dists[$i-$t]){
				     #$max_distance=$clustg_dists[$i-$t];
				 #}#end if
                            #}#end elsif
                            #elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]>2000){
				#$t=-22+$var;
			    #}#end elsif
                            #if($cpt>=6){
				#$cpt=0; 
				#$t=-22+$var;
			    #}#end if
                        #}#end elsif
                }#end for
		#$cpt=0;    
            }#end elsif
	}#end if
  }#end for
}#end for

###############################################################################################################################
for(my $i=0; $i<=scalar(@clustg_ids)-1; $i++){
        if(($clustenzyme_ids[$i] eq $clustg_ids[$i]) && ($positions[$i] == 0)){
            #my @distances;
            #my @max_min;
            #my $max_distance;
            my $cpt=0;
	    #print "POS: $positions[$i]\n";
	    for(my $t=1; $t<=20; $t++){
                #print "$t\n";
		#push(@distances,$clustg_dists[$i-$t]);
		#@max_min = min_and_max(@distances);
		#$max_distance = $max_min[1];
                #print "$max_distance\n";
		if(not defined $clustenzyme_ids[$i-$t]){
		    next;
		}
		if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i-$t]){#INGORE
		    #$t=21;
		}#end if
		elsif($positions[$i-$t]>=0){
		    #print "$positions[$i-$t]\n";
		    if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i-$t]){
			$t=21;
		    }#end if
		    #print "$t\n";
		    if($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]<=3814){
			#print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
			#$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{0}{$domainchecks[$i]}=0;
                        #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
	                $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
                        $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i-$t]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
                        #if($max_distance<$clustg_dists[$i-$t]){
			    #$max_distance=$clustg_dists[$i-$t];
			#}#end if
		    }#end if
			elsif($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]>=3814){
			    #$t=21;
			    $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
			    $t=21;
			}#end elsif
                        elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]<=3814){
			    $cpt=$cpt+1;
                            #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
			    #$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{0}{$domainchecks[$i]}=0;
			    #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                             $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
                             $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i-$t]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
                            #print "$cpt\n";
                            #if($max_distance<$clustg_dists[$i-$t]){
				#$max_distance=$clustg_dists[$i-$t];
			    #}#end if
                        }#end elsif
                        elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]>=3814){
			    #$t=21;
			    $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
			    $t=21;
			}#end elsif
			else{
			}#end else;
                  if($cpt>=10){
			    $cpt=0; 
			    $t=21;
			}#end if
                }#end elsif
		#elsif($positions[$i-$t]<=0){#THIS CONDITIONAL FAILS BECAUSE $positions[$i-$t] is ALWAYS POSITIVE
		    #print "$positions[$i+$t]\n";
		    #if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){
			#$t=21;
		    #}#end if
		    #print "$t\n";
                    #if($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]<2000){
                        #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
			#$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{0}{$domainchecks[$i]}=0;
                        #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                         #$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$gene_order[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
                        #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$gene_order[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0; 
                        #if($max_distance<$clustg_dists[$i+$t]){
			    #$max_distance=$clustg_dists[$i+$t];
			#}#end if
                    #}#end if
                        #elsif($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]>2000){
			    #$t=21;
			#}#end elsif
                        #elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]<2000){
                            #$cpt=$cpt+1;
                            #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                            #$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{0}{$domainchecks[$i]}=0;
			    #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                             #$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$gene_order[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
                             #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$gene_order[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
                            #print "$cpt\n";
                            #if($max_distance<$clustg_dists[$i+$t]){
				#$max_distance=$clustg_dists[$i+$t];
			    #}#end if
                        #}#end elsif
                        #elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]>2000){
			    #$t=21;
			#}#end elsif
                        #else{
			#}#end else;
                        #if($cpt>=6){
			    #$cpt=0; 
			    #$t=21;
			#}#end if
                #}#end elsif
	    }#end for
	    $cpt=0;	    
	    for(my $t=1; $t<=20; $t++){
		#print "$t\n";
                #push(@distances,$clustg_dists[$i-$t]);
                #@max_min = min_and_max(@distances);
                #$max_distance = $max_min[1];
                #print "$max_distance\n";
		if(not defined $clustenzyme_ids[$i+$t]){
		    next;
		}
                if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){#IGNORE
                    #$t=21;
                }#end if
		elsif($positions[$i+$t]<=0){
		    #print "$positions[$i+$t]\n";
		    if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i+$t]){
			$t=21;
		    }#end if
		    #print "$t\n";
                    if($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]<=3814){
                        #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
			#$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{0}{$domainchecks[$i]}=0;
                        #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                         $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
                        $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$t]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;                      
                        #if($max_distance<$clustg_dists[$i+$t]){
			    #$max_distance=$clustg_dists[$i+$t];
			#}#end if
                    }#end if
                        elsif($domainchecks[$i+$t]==1 && $clustg_dists[$i+$t]>=3814){
			    #$t=21;
			    $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
			    $t=21;
			}#end elsif
                        elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]<=3814){
                            $cpt=$cpt+1;
                            #print "$clustenzyme_ids[$i]\t$clustg_ids[$i+$t]\t$positions[$i+$t]\t$clustg_dists[$i+$t]\t$domainchecks[$i+$t]\n";
                            #$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{0}{$domainchecks[$i]}=0;
			    #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}=0;
                             $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
                              $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i+$t]}{$clustg_ids[$i+$t]}{$positions[$i+$t]}{$contigs[$i+$t]}{$starts_5[$i+$t]}{$ends_3[$i+$t]}{$clustg_dists[$i+$t]}{$domainchecks[$i+$t]}{$clustg_functions[$i+$t]}=0;
                            #print "$cpt\n";
                            #if($max_distance<$clustg_dists[$i+$t]){
				#$max_distance=$clustg_dists[$i+$t];
			    #}#end if
                        }#end elsif
                        elsif($domainchecks[$i+$t]==0 && $clustg_dists[$i+$t]>=3814){
			    #$t=21;
			    $clustersn{$clustenzyme_ids[$i]}{$gene_order[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
			    $t=21;
			}#end elsif
                        else{
			}#end else;
                  if($cpt>=10){
			    $cpt=0; 
			    $t=21;
			}#end if
                }#end elsif    
		#elsif($positions[$i+$t]>=0){#THIS CONDITION FAILS BECAUSE $positions[$i+$t] ALWAYS NEGATIVE
		    #print "$t\n";
		    #print "$positions[$i-$t]\n";
		    #if($clustenzyme_ids[$i] !~ $clustenzyme_ids[$i-$t]){
			#$t=21;
		    #}#end if
		    #if($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]<2000){
			#print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
			#$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{0}{$domainchecks[$i]}=0;
                        #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                         #$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$gene_order[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
                         #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$gene_order[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
                        #if($max_distance<$clustg_dists[$i-$t]){
			    #$max_distance=$clustg_dists[$i-$t];
			#}#end if
                    #}#end if
			#elsif($domainchecks[$i-$t]==1 && $clustg_dists[$i-$t]>2000){
			    #$t=21;
			#}#end elsif
			#elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]<2000){
			    #$cpt=$cpt+1;
                            #print "$clustenzyme_ids[$i]\t$clustg_ids[$i-$t]\t$positions[$i-$t]\t$clustg_dists[$i-$t]\t$domainchecks[$i-$t]\n";
                            #$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{0}{$domainchecks[$i]}=0;
			    #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}=0;
                             #$clustersn{$clustenzyme_ids[$i]}{$clustenzyme_ids[$i]}{0}{$contigs[$i]}{$gene_order[$i]}{$starts_5[$i]}{$ends_3[$i]}{0}{$domainchecks[$i]}{$clustg_functions[$i]}=0;
                             #$clustersn{$clustenzyme_ids[$i]}{$clustg_ids[$i-$t]}{$positions[$i-$t]}{$contigs[$i-$t]}{$gene_order[$i-$t]}{$starts_5[$i-$t]}{$ends_3[$i-$t]}{$clustg_dists[$i-$t]}{$domainchecks[$i-$t]}{$clustg_functions[$i-$t]}=0;
                            #print "$cpt\n";
			    #if($max_distance<$clustg_dists[$i-$t]){
				#$max_distance=$clustg_dists[$i-$t];
			    #}#end if
			#}#end elsif
			#elsif($domainchecks[$i-$t]==0 && $clustg_dists[$i-$t]>2000){
			    #$t=21;
			#}#end elsif
		        #else{
			#}#end else;
		        #if($cpt>=6){
			    #$cpt=0; 
			    #$t=21;
			#}#end if
	       #}#end elsif 
	    }#end for
            $cpt=0;		
	}#end if
}#end for
############################################################################################################################################

#for(my $t=-6; $t>-20; --$t){
    #print "$t\n";	
#}

#foreach my $c_e(sort keys %clustersn){
    #foreach my $c_g(sort keys %{$clustersn{$c_e}}){
	#foreach my $pos(sort keys %{$clustersn{$c_e}{$c_g}}){
	    #foreach my $dis(sort keys %{$clustersn{$c_e}{$c_g}{$pos}}){
		#foreach my $dom(sort keys %{$clustersn{$c_e}{$c_g}{$pos}{$dis}}){
		    #print "$c_e\t$c_g\t$pos\t$dis\t$dom\n";
		#}#end foreach5
	   #}#end foreach4
	#}#end foreach3
    #}#end foreach2
#}#end foreach1

foreach my $c_e(sort keys %clustersn){
    foreach my $c_g(sort keys %{$clustersn{$c_e}}){
        foreach my $pos(sort keys %{$clustersn{$c_e}{$c_g}}){
	    foreach my $contigs(sort keys %{$clustersn{$c_e}{$c_g}{$pos}}){
		foreach my $gene_order(sort keys %{$clustersn{$c_e}{$c_g}{$pos}{$contigs}}){
		    foreach my $start_5(sort keys %{$clustersn{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}}){
			foreach my $end_3(sort keys %{$clustersn{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}{$start_5}}){
                            foreach my $dis(sort keys %{$clustersn{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}{$start_5}{$end_3}}){
				foreach my $dom(sort keys %{$clustersn{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}{$start_5}{$end_3}{$dis}}){
				    foreach my $g_function(sort keys %{$clustersn{$c_e}{$c_g}{$pos}{$contigs}{$gene_order}{$start_5}{$end_3}{$dis}{$dom}}){
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
						
    foreach my $c_e(sort keys %clustersn){
	foreach my $gene_order(sort {$a <=> $b} keys %{$clustersn{$c_e}}){
	    foreach my $g_id(sort keys %{$clustersn{$c_e}{$gene_order}}){
		foreach my $pos(sort keys %{$clustersn{$c_e}{$gene_order}{$g_id}}){
		    foreach my $chromos(sort keys %{$clustersn{$c_e}{$gene_order}{$g_id}{$pos}}){
			foreach my $start_5(sort keys %{$clustersn{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}}){
			    foreach my $end_3(sort keys %{$clustersn{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}}){
				foreach my $g_distance(sort keys %{$clustersn{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}}){
				    foreach my $domain_score(sort keys %{$clustersn{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}{$g_distance}}){
					foreach my $gene_function(sort keys %{$clustersn{$c_e}{$gene_order}{$g_id}{$pos}{$chromos}{$start_5}{$end_3}{$g_distance}{$domain_score}}){
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
exit 0;

########################################################################
#
#subroutine to open files, filenames provided as command line arguments
#
########################################################################
sub open_file {
  my($filename, $header_list, $sequence_list) = @_;
  my $line; # read one-line at time from <INFILE>
  my $counter = -1; # iterate through all headers and sequences in <INFILE>

  # open file
  open(INFILE, $filename) or die "Can't open $filename";
  while ($line = <INFILE>) {
    chomp $line;                         # Remove extra end line \n from EOF

    # If loop to skip to the next line if the first line in file is blank
    if($line =~ /^\s*$/) {
      next;
    } #end if
    elsif ($line =~ /^>/) {              # Line starts with a ">"
      $counter++;                        # Increament counter to go to next line
      $$header_list[$counter]=($line);   # Save header line
      $$sequence_list[$counter] = '';    # Start with an empty sequence
    } #end elsif
    else {
      $$sequence_list[$counter]=$$sequence_list[$counter] . $line;     # Add line to end of sequence
      $$sequence_list[$counter] =~ s/(\s|\n)//g;   # Remove all white space and end line (\n) from sequence
    } #end else
  } #end while
} #end sub open_file
########################################################################

########################################################################
#
#subroutine to select minimum and maximum gene_distances
#
########################################################################
sub min_and_max {
  my (@numbers);
  @numbers = @_;
  my ($min, $max);
  $min = $numbers[0];
  $max = $numbers[0];

  foreach my $i (@numbers) {
    if($i > $max) {
      $max = $i;
    } # end if($i > $max)
    elsif($i < $min) {
      $min = $i;
    } # end elsif($i < $min)
  } # end foreach my $i (@numbers)

  return ($min, $max);
}#end subroutine min_and_max
########################################################################
