#!/usr/bin/env perl -w

# Author: Fayaz Seifuddin, Nora Khaldi
# Date Created: 06/12/2007
# Last Date Modified: 06/12/2007
# File:testscript4.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: Read in CLUSTERINFOOUT file, perform HMMsearches for domains most commonly found in A.fumigatus clusters
# parse HMMsearches , STEP4
# Usage: perl testscript4.pl GENEDISTOUT $tempdir

use strict;

# Variable declaration
my $filename = $ARGV[0];  # CLUSTERINFOOUT file
my $tempdir = $ARGV[1]; #$tempdir to store HMMoutput
my $base_url = 'data/hmmsgeneclusters';

# Rum HMM searches
print "-- (01 / 31) Searching AA_permease.hmm...";
system "hmmsearch --cut_tc $base_url/AA_permease.hmm $filename > $tempdir/hmm_AA_permease.output";
print "done!\n";

print "-- (02 / 31) Searching ABC_membrane.hmm...";
system "hmmsearch --cut_tc $base_url/ABC_membrane.hmm $filename > $tempdir/ABC_membrane.output";
print "done!\n";

print "-- (03 / 31) Searching ABC_tran.hmm...";
system "hmmsearch --cut_tc $base_url/ABC_tran.hmm $filename > $tempdir/ABC_tran.output";
print "done!\n";

print "-- (04 / 31) Searching ADH_zinc_N.hmm...";
system "hmmsearch --cut_tc $base_url/ADH_zinc_N.hmm $filename > $tempdir/ADH_zinc_N.output";
print "done!\n";

print "-- (05 / 31) Searching Acetyltransf_1.hmm...";
system "hmmsearch --cut_tc $base_url/Acetyltransf_1.hmm $filename > $tempdir/Acetyltransf_1.output";
print "done!\n";

print "-- (06 / 31) Searching FAD_binding_3.hmm...";
system "hmmsearch --cut_tc $base_url/FAD_binding_3.hmm $filename > $tempdir/FAD_binding_3.output";
print "done!\n";

print "-- (07 / 31) Searching FAD_binding_4.hmm...";
system "hmmsearch --cut_tc $base_url/FAD_binding_4.hmm $filename > $tempdir/FAD_binding_4.output";
print "done!\n";

print "-- (08 / 31) Searching Fungal_trans.hmm...";
system "hmmsearch --cut_tc $base_url/Fungal_trans.hmm $filename > $tempdir/Fungal_trans.output";
print "done!\n";

print "-- (09 / 31) Searching MFS_1.hmm...";
system "hmmsearch --cut_tc $base_url/MFS_1.hmm $filename > $tempdir/MFS_1.output";
print "done!\n";

print "-- (10 / 31) Searching Methyltransf_2.hmm...";
system "hmmsearch --cut_tc $base_url/Methyltransf_2.hmm $filename > $tempdir/Methyltransf_2.output";
print "done!\n";

print "-- (11 / 31) Searching Pyr_redox.hmm...";
system "hmmsearch --cut_tc $base_url/Pyr_redox.hmm $filename > $tempdir/Pyr_redox.output";
print "done!\n";

print "-- (12 / 31) Searching Sugar_tr.hmm...";
system "hmmsearch --cut_tc $base_url/Sugar_tr.hmm $filename > $tempdir/Sugar_tr.output";
print "done!\n";

print "-- (13 / 31) Searching Transferase.hmm...";
system "hmmsearch --cut_tc $base_url/Transferase.hmm $filename > $tempdir/Transferase.output";
print "done!\n";

print "-- (14 / 31) Searching Zn_clus.hmm...";
system "hmmsearch --cut_tc $base_url/Zn_clus.hmm $filename > $tempdir/Zn_clus.output";
print "done!\n";

print "-- (15 / 31) Searching adh_short.hmm...";
system "hmmsearch --cut_tc $base_url/adh_short.hmm $filename > $tempdir/adh_short.output";
print "done!\n";

print "-- (16 / 31) Searching p450.hmm...";
system "hmmsearch --cut_tc $base_url/p450.hmm $filename > $tempdir/p450.output";
print "done!\n";

#Adding new domains#
print "-- (17 / 31) Searching GST_C.hmm...";
system "hmmsearch --cut_tc $base_url/GST_C.hmm $filename > $tempdir/GST_C.output";
print "done!\n";

print "-- (18 / 31) Searching GST_N.hmm...";
system "hmmsearch --cut_tc $base_url/GST_N.hmm $filename > $tempdir/GST_N.output";
print "done!\n";

print "-- (19 / 31) Searching Peptidase_M19.hmm...";
system "hmmsearch --cut_tc $base_url/Peptidase_M19.hmm $filename > $tempdir/Peptidase_M19.output";
print "done!\n";

print "-- (20 / 31) Searching 2OG-FeII_Oxy.hmm...";
system "hmmsearch --cut_tc $base_url/2OG-FeII_Oxy.hmm $filename > $tempdir/2OG-FeII_Oxy.output";
print "done!\n";

print "-- (21 / 31) Searching PhyH.hmm...";
system "hmmsearch --cut_tc $base_url/PhyH.hmm $filename > $tempdir/PhyH.output";
print "done!\n";

print "-- (22 / 31) Searching Pyr_redox_2.hmm...";
system "hmmsearch --cut_tc $base_url/Pyr_redox_2.hmm $filename > $tempdir/Pyr_redox_2.output";
print "done!\n";

print "-- (23 / 31) Searching Aminotran_1_2.hmm...";
system "hmmsearch --cut_tc $base_url/Aminotran_1_2.hmm $filename > $tempdir/Aminotran_1_2.output";
print "done!\n";

print "-- (24 / 31) Searching Lactamase_B.hmm...";
system "hmmsearch --cut_tc $base_url/Lactamase_B.hmm $filename > $tempdir/Lactamase_B.output";
print "done!\n";

print "-- (25 / 31) Searching Peptidase_M24.hmm...";
system "hmmsearch --cut_tc $base_url/Peptidase_M24.hmm $filename > $tempdir/Peptidase_M24.output";
print "done!\n";

# Removed --cut_tc : use Pfam TC trusted threshold cutoffs
print "-- (26 / 31) Searching GTP_cyclohydroI.hmm...";
system "hmmsearch $base_url/GTP_cyclohydroI.hmm $filename > $tempdir/GTP_cyclohydroI.output";
print "done!\n";

print "-- (27 / 31) Searching methyl_EasF.hmm...";
system "hmmsearch $base_url/methyl_EasF.hmm $filename > $tempdir/methyl_EasF.output";
print "done!\n";


print "-- (28 / 31) Searching Epimerase.hmm...";
system "hmmsearch $base_url/Epimerase.hmm $filename > $tempdir/Epimerase.output";
print "done!\n";

print "-- (29 / 31) Searching ergot_EASG.hmm...";
system "hmmsearch $base_url/ergot_EASG.hmm $filename > $tempdir/ergot_EASG.output";
print "done!\n";

print "-- (30 / 31) Searching Oxidored_FMN.hmm...";
system "hmmsearch $base_url/Oxidored_FMN.hmm $filename > $tempdir/Oxidored_FMN.output";
print "done!\n";

print "-- (31 / 31) Searching tRNA-synt_2.hmm...";
system "hmmsearch $base_url/tRNA-synt_2.hmm $filename > $tempdir/tRNA-synt_2.output";
print "done!\n";

exit;
