#!/usr/bin/env perl -w

use warnings;
use strict;
use File::Basename;
my $dirname = dirname(__FILE__);

# Variable Declaration
my $infile = $ARGV[0]; #protein input
my $infile1 = $ARGV[1]; #coordinates file
my $tempdir = $ARGV[2];


my $baseurl = 'pl/';

my $create_input = $baseurl.'create_genome_input.pl';
my $hmm_search_nrps_pks_dmat = $baseurl.'HMM_search_nrps_pks_dmat.pl';
my $nrps_pks_dmat_content = $baseurl.'determine_nrps_pks_dmat.pl';
my $make_clusters = $baseurl.'make_SM_clusters.pl';
my $hmm_search_clusters = $baseurl.'HMM_search_clusters.pl';
my $clusters_domain_content = $baseurl.'clusters_domain_content.pl';
my $sort_clusters_by_gene_order = $baseurl.'sort_clusters_by_gene_order.pl';
my $first_delimited_clusters = $baseurl.'create_clusters_gdistance_domaincontent_version2_08.07.2008.pl'; #changed 08.07.08 to version2 script
my $final_clusters = $baseurl.'delimit_clusters_further_version3_08.07.2008.pl'; #changed 08.07.08 to version3 script
my $non_duplicate_clusters = $baseurl.'SMURF_duplicate_blast.pl';#added new script to remove duplicate clusters 08.12.08
my $nrps_pks_dmat_content_attachment = $baseurl.'get_nrps_pks_dmat_file_user.pl'; #changed 08.12.08 to original script
my $htab = $baseurl.'htab.pl';

my $tmp;

# run the program on the input file and save the output
system "mkdir -p $tempdir";
system "chmod 777 $tempdir";

print "(01 / 11) Executing $create_input\n";
system "$create_input $infile $infile1 > $tempdir/OUT.input";

print "(02 / 11) Executing $hmm_search_nrps_pks_dmat\n";
system "$hmm_search_nrps_pks_dmat $tempdir/OUT.input $tempdir";

# Parse HMMsearches outputs
#$tmp = `ls $tempdir/*.out | $htab -f`;
$tmp = `python $dirname/jcvi_htab.py $tempdir '.out'`;

print "(03 / 11) Executing $nrps_pks_dmat_content\n";
system "$nrps_pks_dmat_content $tempdir/OUT.input $tempdir > $tempdir/DOMAINCONTENTOUT.input";

print "(04 / 11) Executing $nrps_pks_dmat_content_attachment\n";
system "$nrps_pks_dmat_content_attachment $tempdir/OUT.input $tempdir > $tempdir/NRPS_PKS_DMAT_ATTACHMENT.txt";

print "(05 / 11) Executing $make_clusters\n";
system "$make_clusters $tempdir/OUT.input $tempdir/DOMAINCONTENTOUT.input > $tempdir/GENEDISTOUT.input";

print "(06 / 11) Executing $hmm_search_clusters\n";
system "$hmm_search_clusters $tempdir/GENEDISTOUT.input $tempdir";

# Parse HMMsearches outputs
#$tmp = `ls $tempdir/*.output | $htab -f`;
$tmp = `python $dirname/jcvi_htab.py $tempdir '.output'`; #this is currently failing though, not sure why yet

print "(07 / 11) Executing $clusters_domain_content\n";
system "$clusters_domain_content $tempdir/GENEDISTOUT.input $tempdir > $tempdir/CLUSTER-DOMAIN-CONTENT-OUT.input";

print "(08 / 11) Executing $sort_clusters_by_gene_order\n";
system "$sort_clusters_by_gene_order $tempdir/CLUSTER-DOMAIN-CONTENT-OUT.input > $tempdir/NEW-CLUSTER-DOMAIN-CONTENT-OUT.input";

print "(09 / 11) Executing $first_delimited_clusters\n";
system "$first_delimited_clusters $tempdir/DOMAINCONTENTOUT.input $tempdir/NEW-CLUSTER-DOMAIN-CONTENT-OUT.input > $tempdir/FIRST_DELIMITED_CLUSTERS.output";

print "(10 / 11) Executing $final_clusters\n";
system "$final_clusters $tempdir/FIRST_DELIMITED_CLUSTERS.output > $tempdir/FINAL_CLUSTERS.output";

print "(11 / 11) Executing $non_duplicate_clusters\n";
system "$non_duplicate_clusters $tempdir/FINAL_CLUSTERS.output > $tempdir/NON_DUPLICATE_FINAL_CLUSTERS.txt";

system "chmod -R 777 $tempdir";