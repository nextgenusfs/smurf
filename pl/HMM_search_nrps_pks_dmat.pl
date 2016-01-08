#!/usr/bin/env perl -w

# Author: Fayaz Seifuddin, Nora Khaldi
# Date Created: 05/22/2007
# Last Date Modified: 05/23/2007
# File:testscript2.pl
# Sources of Code: Originally created by Fayaz Seifuddin, Nora Khaldi
# Program Description: Read in inputFASTA file, perform hmmsearches, Pasre hmmsearches,(NRPS&PKS domains),STEP2
# Usage: perl testscript2.pl inputFASTA file $tempdir

use strict;

# Variable declaration
my $filename = $ARGV[0];  # Input fasta file
my $tempdir = $ARGV[1]; # Temporary directory to store HMMoutput
my $base_url = 'data/hmmsPKSNRPS';

# Run HMM searches
print "-- (01 / 13) Searching AMP-binding.hmm...";
system "hmmsearch --cut_tc $base_url/AMP-binding.hmm $filename > $tempdir/hmm_AMP.out";
print "done!\n";

print "-- (02 / 13) Searching Acyl_transf_1.hmm...";
system "hmmsearch --cut_tc $base_url/Acyl_transf_1.hmm  $filename > $tempdir/hmm_Acyl_transf.out";
print "done!\n";

print "-- (03 / 13) Searching Condensation.hmm...";
system "hmmsearch --cut_tc $base_url/Condensation.hmm $filename > $tempdir/hmm_Condensation.out";
print "done!\n";

print "-- (04 / 13) Searching Ketoacyl-synt_C.hmm...";
system "hmmsearch --cut_tc $base_url/Ketoacyl-synt_C.hmm $filename > $tempdir/hmm_Ketoacyl-synt_C.out";
print "done!\n";

print "-- (05 / 13) Searching PP-binding.hmm...";
system "hmmsearch --cut_tc $base_url/PP-binding.hmm $filename > $tempdir/hmm_PP-binding.out";
print "done!\n";

print "-- (06 / 13) Searching ketoacyl-synt.hmm...";
system "hmmsearch --cut_tc $base_url/ketoacyl-synt.hmm $filename > $tempdir/hmm_ketoacyl-synt.out";
print "done!\n";

# Adding 2 new domains for NRPS-like enzymes
print "-- (07 / 13) Searching NAD_binding_4.hmm...";
system "hmmsearch --cut_tc $base_url/NAD_binding_4.hmm $filename > $tempdir/hmm_NAD_binding_4.out";
print "done!\n";

print "-- (08 / 13) Searching adh_short.hmm...";
system "hmmsearch --cut_tc $base_url/adh_short.hmm $filename > $tempdir/hmm_adh_short.out";
print "done!\n";

# Adding the DMAT search domains
print "-- (09 / 13) Searching DMAT.HMM...";
system "hmmsearch $base_url/DMAT.HMM $filename > $tempdir/hmm_DMAT.out";
print "done!\n";

print "-- (10 / 13) Searching arom_pren_DMATS.hmm...";
system "hmmsearch $base_url/arom_pren_DMATS.hmm $filename > $tempdir/hmm_arom_pren_DMATS.out";
print "done!\n";

print "-- (11 / 13) Searching trp_dimet_allyl.hmm...";
system "hmmsearch $base_url/trp_dimet_allyl.hmm $filename > $tempdir/hmm_trp_dimet_allyl.out";
print "done!\n";

# Adding domain to eliminate false positive rate for C-terminal domain of L-aminoadipate-semialdehyde dehydrogenase alpha subunit
print "-- (12 / 13) Searching Lys2.hmm...";
system "hmmsearch $base_url/Lys2.hmm $filename > $tempdir/hmm_Lys2.out";
print "done!\n";

print "-- (13 / 13) Searching alpha_am_amid.hmm...";
system "hmmsearch $base_url/alpha_am_amid.hmm $filename > $tempdir/hmm_alpha_am_amid.out";
print "done!\n";

exit;