#!/usr/local/bin/perl -w 

use strict;

# Date 08/08/08
# This script takes the SMURF output file. Gets rid of duplicates and merges others
#To run it ./SMURF_duplicate_blast.pl SMURF_Final_Cluster_Output.tab
#The output of this script is a tab file named Final_Output_No_Duplicates_SMURF.tab

##################################################  Variables ###################
my %hash;
my %hash_all;
my %array_backbones;
my $i=0;
my $cpt=0;
#my $OUTFILE;
my $backbone;


################################################### Reading the input SMURF file ####################################
local $/= "\n#";

open(IN,$ARGV[0]);# File with SMURF cluster output
	while (<IN>) {
	$cpt=0;
	my @line=split(/\n/,$_); 
		if($cpt==0 ){my @characters=split(/\t/,$line[0]); $backbone=$characters[0];$cpt=1;}
		for ($i=0;$i<$#line;$i++){
			my @characters=split(/\t/,$line[$i]);
			$hash{$backbone}{$characters[1]}=$characters[2];
			$hash_all{$backbone}{$characters[1]}=$line[$i];	
		}
		$array_backbones{$.}=$backbone;

	}
close IN;

################################################### Testing the existance of duplicates and larger clusters  ####################################

foreach my $k(sort keys %array_backbones){
	foreach my $k1(sort keys %hash){
		foreach my $k2 (sort keys %{$hash{$k1}}){
			if($k2 eq $array_backbones{$k} && $k1 ne $array_backbones{$k}){# I check if the backbone is in teh cluster without being the main gene of the cluster.
				if (keys(%{$hash{$k1}}) eq keys(%{$hash{$array_backbones{$k}}})) {# If the two clusters are identical I only keep one.
					&deleting(\%hash,$array_backbones{$k});	
				}
				else{# If the cluters are different I will take the biggest and add to it anything in the other.
					if (&one_in_other(\%hash,$array_backbones{$k},$k1)==0){
						if(keys(%{$hash{$k1}}) > keys(%{$hash{$array_backbones{$k}}})){
							&deleting(\%hash,$array_backbones{$k});		
						}
						elsif(keys(%{$hash{$k1}}) < keys(%{$hash{$array_backbones{$k}}})){
							&deleting(\%hash,$k1);
						}	
					}
				}
			}
		}
	}
}

################################################### Printing out the SMURF clusters without the duplicates ####################################


$i=0;
#open OUTFILE,"> Final_Output_No_Duplicates_SMURF.tab" or die "can't open output file Final_Output_No_Duplicates_SMURF.tab\n";
foreach my $k1(sort keys %hash){
	if (keys(%{$hash{$k1}})>1){# this is because I can't get rid of empty elements in the hash. I need to come back to this!!
		$i++;
        	#print OUTFILE "Cluster:$i\n";
		print "Cluster:$i\n";
        	#print OUTFILE "Backbone_gene_id Gene_id Gene_positions  Chromosome-Contig       Gene_order      5'end   3'end   Gene_distance   Domain_score    Annotated_gene_function\n"; 
		print "Backbone_gene_id\tGene_id\tGene_positions\tChromosome-Contig\tGene_order\t5'end\t3'end\tGene_distance\tDomain_score\tAnnotated_gene_function\n";
		foreach my $k2 (sort {$hash{$k1}{$a}<=>$hash{$k1}{$b}} keys %{$hash{$k1}}){
        		#print OUTFILE "$hash_all{$k1}{$k2}\n";
			print "$hash_all{$k1}{$k2}\n";
		}
		#print OUTFILE "\n";
		print "\n";
	}
}
#close OUTFILE;

#print"The number of clusters in '$ARGV[0]' is : $i\n";

################################################### Subroutines  ####################################


## 1st subroutine ##

sub one_in_other{
my $f=$_[1];# this is the bckbone that is in a cluter=$array_backbones{$k}
my $f1=$_[2];# this is the backbone that is teh main backbone=$k1
my %hash=%{$_[0]};
my $cpt=1;
my $nono;
my $nono1;

if (keys(%{$hash{$f1}}) > keys(%{$hash{$f}})){

	foreach my $k3 (sort keys %{$hash{$f1}}){
		if(exists $hash{$k3}{$f1}){
        		$nono.=$hash{$k3}{$f1};
		}
        }

	foreach my $k2 (sort keys %{$hash{$f}}){
		if($hash{$f}{$k2} != $nono){$cpt=0; }
	}
}

else{
        foreach  my $k4 (sort keys %{$hash{$f}}){
                if(exists $hash{$k4}{$f}){
			$nono1.=$hash{$k4}{$f};
		}
        }

        foreach  my $k5 (sort keys %{$hash{$f1}}){
                if($hash{$f1}{$k5} != $nono1){$cpt=0; }
        }

}
return $cpt;
}


## 2nd subroutine ##

sub deleting{
my $f=$_[1];
my %hash=%{$_[0]};
foreach my $k2 (sort keys %{$hash{$f}}){
	delete $hash{$f}{$k2};
}
return %hash;
}



exit 0;
