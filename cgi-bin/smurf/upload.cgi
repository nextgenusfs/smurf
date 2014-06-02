#!/usr/local/bin/perl -w
# Purpose: Perl-CGI web interface for testing
# A simple web interface.
# File: upload.cgi
# Usage: http://natalief-lx:8080/cgi-bin/upload.cgi
# Author: Fayaz Seifuddin
# Sources of Code: Orginally created by Fayaz Seifuddin
# Date Created: 08/23/07
# Last Date Modified: 08/23/07
# cookie-get.cgi - fetch the value of a cookie

$ENV{SYBASE}="/usr/local/packages/sybase";

use warnings;
use strict;
use CGI;
use CGI qw(:standard);
use CGI::Carp  qw/fatalsToBrowser/;
use DBI;
use HTML::Template;
use CGI::Session;
use File::Temp qw/ tempfile tempdir /;
use MIME::Lite;

open STDOUT;
open STDERR;
open STDIN;

# Variable Declaration
my $tmpfile;
my $infile;
my $infile1;
my $infile2;
my $outfile;

my $redir = '/smurf/submit.php';
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

my $q = new CGI;

print $q->header(-type=>'text/html');

my $email = $q->cookie('smurfs');

my $proteinFASTAfile =  $q->param('FASTA');
my $proteinFASTAsequence = $q->param('FASTAseq');
my $coordinatesfile = $q->param('coords');
my $coordinatesequence = $q->param('coordinates');
my $uploadgenome = '/opt/www/smurf/data/outfiles/OUT_'.$q->param('genomes').'.input';

my $proteinFASTAfileHandle = $q->upload('FASTA');
my $coordinatesfileHandle = $q->upload('coords');

if(defined $proteinFASTAfileHandle) {
  # get name of uploaded proteinFASTAfileHandle in local tmp directory
  $tmpfile = tmpFileName($proteinFASTAfileHandle);
  close $proteinFASTAfileHandle;
  
  $infile = "$tmpfile.in";   # name of converted input file
  convert_to_unix($infile, $tmpfile);
} # end if(defined $proteinFASTAfileHandle)

if(defined $coordinatesfileHandle) {
  # get name of uploaded coordinatesFASTAfileHandle in local tmp directory
  $tmpfile = tmpFileName($coordinatesfileHandle);
  close $coordinatesfileHandle;

  $infile1 = "$tmpfile.in";   # name of converted input file
  convert_to_unix($infile1, $tmpfile);
} # end if(defined $coordinatesfileHandle)

if($proteinFASTAsequence) {
  #print "$proteinFASTAsequence\n";

  # get hand input of proteinFASTAsequence in local tmp directory
  $tmpfile= "/var/tmp/prot$$";

  open(TEMP_FILE, ">$tmpfile");
  print TEMP_FILE "$proteinFASTAsequence";
  close(TEMP_FILE);

  $infile = "$tmpfile.in";   # name of converted input file
  convert_to_unix($infile, $tmpfile);
} # end if($proteinFASTAsequence)

if($coordinatesequence) {
  #print "$coordinatesequence\n";
  
  # get hand input of coordinatesequence in local tmp directory
  $tmpfile= "/var/tmp/coord$$";

  open(TEMP_FILE, ">$tmpfile");
  print TEMP_FILE "$coordinatesequence";
  close(TEMP_FILE);

  $infile1 = "$tmpfile.in";   # name of converted input file
  convert_to_unix($infile1, $tmpfile);
} # end if($coordinatesequence)

if($uploadgenome) {
  $infile2 = $uploadgenome;
} # end if($uploadgenome)

# Print the proteinFASTAfile on the screen
#print p, "<b>proteinFASTAfile: </b>";
#print p, "<b>Please page down for more results...</b>";
#print p;
#print_file($infile);
#print hr;

# Print the coordinatesfile on the screen
#print p, "<b>coordinatesfile: </b>";
#print p, "<b>Please page down for more results...</b>";
#print p;
#print_file($infile1);
#print hr;

# Print the inputfile on the screen
#print p, "<b>inputfile: </b>";
#print p, "<b>Please page down for more results...</b>";
#print p;
#print_file($infile2);
#print hr;

my $pid = fork;

if($pid) {
  #print_redirect($redir);
} #end if
else{
  &runSMURF($infile,$infile1,$infile2,$email);
 
  # remove input and output files (the tmpfile is automatically deleted)
  unlink $infile,$infile1,$infile2;
} #end else

print_redirect($redir);
exit;

################################################################################
# Subroutine to runSMURF
################################################################################
sub runSMURF{
  my ($infile, $infile1, $infile2, $email) = @_;
  my $tempdir = tempdir("/opt/www/smurf/tmp/XXXX", CLEANUP=>1);
  my $base_url = '/opt/www/smurf/data/logs';
  my $tmp;

  # run the program on the input file and save the output
  print "This may take several minutes. The results will be sent to you by email. <br>\n";
  print "Questions: <a href=\"mailto:smurf\@jcvi.org\">smurf\@jcvi.org</a><br><br>\n";

  if(!($infile2 eq "/opt/www/smurf/data/outfiles/OUT_.input")) {
    system "chmod 777 $tempdir";
    
    print "(1 / 10) Executing $hmm_search_nrps_pks_dmat <br>\n";
    system "$hmm_search_nrps_pks_dmat $infile2 $tempdir";

    # Parse HMMsearches outputs
    $tmp = `ls $tempdir/*.out | $htab -f`;

    print "(2 / 10) Executing $nrps_pks_dmat_content <br>\n";
    system "$nrps_pks_dmat_content $infile2 $tempdir > $tempdir/DOMAINCONTENTOUT.input";

    print "(3 / 10) Executing $nrps_pks_dmat_content_attachment <br>\n";
    system "$nrps_pks_dmat_content_attachment $infile2 $tempdir > $tempdir/NRPS_PKS_DMAT_ATTACHMENT.txt";

    print "(4 / 10) Executing $make_clusters <br>\n";
    system "$make_clusters $infile2 $tempdir/DOMAINCONTENTOUT.input > $tempdir/GENEDISTOUT.input";

    print "(5 / 10) Executing $hmm_search_clusters <br>\n";
    system "$hmm_search_clusters $tempdir/GENEDISTOUT.input $tempdir";

    # Parse HMMsearches outputs
    $tmp = `ls $tempdir/*.output | $htab -f`;

    print "(6 / 10) Executing $clusters_domain_content <br>\n";
    system "$clusters_domain_content $tempdir/GENEDISTOUT.input $tempdir > $tempdir/CLUSTER-DOMAIN-CONTENT-OUT.input";

    print "(7 / 10) Executing $sort_clusters_by_gene_order <br>\n";
    system "$sort_clusters_by_gene_order $tempdir/CLUSTER-DOMAIN-CONTENT-OUT.input > $tempdir/NEW-CLUSTER-DOMAIN-CONTENT-OUT.input";

    print "(8 / 10) Executing $first_delimited_clusters <br>\n";
    system "$first_delimited_clusters $tempdir/DOMAINCONTENTOUT.input $tempdir/NEW-CLUSTER-DOMAIN-CONTENT-OUT.input > $tempdir/FIRST_DELIMITED_CLUSTERS.output";

    print "(9 / 10) Executing $final_clusters <br>\n";
    system "$final_clusters $tempdir/FIRST_DELIMITED_CLUSTERS.output > $tempdir/FINAL_CLUSTERS.output";

    print "(10 / 10) Executing $non_duplicate_clusters <br>\n";
    system "$non_duplicate_clusters $tempdir/FINAL_CLUSTERS.output > $tempdir/NON_DUPLICATE_FINAL_CLUSTERS.txt";

    system "chmod -R 777 $tempdir";
  } #end if
  else{
    # run the program on the input file and save the output
    system "chmod 777 $tempdir";

    print "(01 / 11) Executing $create_input<br>";
    system "$create_input $infile $infile1 > $tempdir/OUT.input";

    print "(02 / 11) Executing $hmm_search_nrps_pks_dmat<br>";
    system "$hmm_search_nrps_pks_dmat $tempdir/OUT.input $tempdir";
    
    # Parse HMMsearches outputs
    $tmp = `ls $tempdir/*.out | $htab -f`;
    
    print "(03 / 11) Executing $nrps_pks_dmat_content<br>";
    system "$nrps_pks_dmat_content $tempdir/OUT.input $tempdir > $tempdir/DOMAINCONTENTOUT.input";
    
    print "(04 / 11) Executing $nrps_pks_dmat_content_attachment<br>";
    system "$nrps_pks_dmat_content_attachment $tempdir/OUT.input $tempdir > $tempdir/NRPS_PKS_DMAT_ATTACHMENT.txt";
    
    print "(05 / 11) Executing $make_clusters";
    system "$make_clusters $tempdir/OUT.input $tempdir/DOMAINCONTENTOUT.input > $tempdir/GENEDISTOUT.input";
    
    print "(06 / 11) Executing $hmm_search_clusters<br>";
    system "$hmm_search_clusters $tempdir/GENEDISTOUT.input $tempdir";    

    # Parse HMMsearches outputs
    $tmp = `ls $tempdir/*.output | $htab -f`;

    print "(07 / 11) Executing $clusters_domain_content<br>";
    system "$clusters_domain_content $tempdir/GENEDISTOUT.input $tempdir > $tempdir/CLUSTER-DOMAIN-CONTENT-OUT.input";

    print "(08 / 11) Executing $sort_clusters_by_gene_order<br>";
    system "$sort_clusters_by_gene_order $tempdir/CLUSTER-DOMAIN-CONTENT-OUT.input > $tempdir/NEW-CLUSTER-DOMAIN-CONTENT-OUT.input";
    
    print "(09 / 11) Executing $first_delimited_clusters<br>";
    system "$first_delimited_clusters $tempdir/DOMAINCONTENTOUT.input $tempdir/NEW-CLUSTER-DOMAIN-CONTENT-OUT.input > $tempdir/FIRST_DELIMITED_CLUSTERS.output";

    print "(10 / 11) Executing $final_clusters<br>";
    system "$final_clusters $tempdir/FIRST_DELIMITED_CLUSTERS.output > $tempdir/FINAL_CLUSTERS.output";
    
    print "(11 / 11) Executing $non_duplicate_clusters <br>\n";
    system "$non_duplicate_clusters $tempdir/FINAL_CLUSTERS.output > $tempdir/NON_DUPLICATE_FINAL_CLUSTERS.txt";

    system "chmod -R 777 $tempdir";
  }#end else

  my $from_address = 'smurf@jcvi.org';
  my $to_address = $email;
  my $bcc_address = 'smurf@jcvi.org';
  my $mail_host = 'mailhost.tigr.org';

  my $subject = 'SMURF Output';
  my $message_body = 'Dear researcher,

    Attached you will find your Secondary Metabolite clusters & Backbone genes outputs from SMURF.
    
    If you have any questions or comments, please contact us at smurf@jcvi.org         
    
    Privacy note: We will not share your confidential data with anyone. The uploaded files were kept briefly on our secure server and permanently deleted.

    Sincerely,

    SMURF team';

  my $file_path = $tempdir.'/NON_DUPLICATE_FINAL_CLUSTERS.txt';
  my $file_name = 'Secondary-Metabolite-Clusters.txt';
  my $file_path2 = $tempdir.'/NRPS_PKS_DMAT_ATTACHMENT.txt';
  my $file_name2 = 'Backbone-genes.txt';

  my $msg = MIME::Lite->new(
    From=>$from_address,
    To=>$to_address,
    Bcc=>$bcc_address,  
    Subject=>$subject,
    Type=>'multipart/mixed'
  )or die "Error creating multipart container: $!\n";

  $msg->attach(
    Type=>'TEXT',
    Data=>$message_body
  )or die "Error adding the text message part: $!\n";

  $msg->attach(
    Type=>'text/plain',
    Path=>$file_path,
    Filename=>$file_name,
    Disposition=>'attachment'
  )or die "Error attaching $file_name: $!\n";

  $msg->attach(
    Type=>'text/plain',
    Path=>$file_path2,
    Filename=>$file_name2,
    Disposition=>'attachment'
  )or die "Error attaching $file_name2: $!\n";

  MIME::Lite->send('smtp',$mail_host);
  $msg->send() or die('Error sending message: $!');

  close STDOUT;
  close STDERR;
  close STDIN;
}# end runSMURF
################################################################################

################################################################################
# Subroutine to redirect web page to any URL after the user hits the submit 
# button
################################################################################
sub print_redirect {
  my ($url) = @_;

  # get the redirector template
  my $tmpl = HTML::Template->new( filename => '../../htdocs/smurf/templates/redirect.tmpl',
                                  die_on_bad_params => 1
				);

  $tmpl->param( REDIR => $url );
  print $tmpl->output;
} # end print_redirect
################################################################################

################################################################################
# Subroutine to print PREFORMATTED OUTPUT on the screen
################################################################################
sub print_file {
  my ($file) = @_;

  if (open(OUTFILE, "$file")) {
    my @output = <OUTFILE>;
    close OUTFILE;

    print "<PRE>";              # preformatted output
    
    foreach my $line (@output) {
      # avoid printing special HTML characters
      $line =~ s/&/&amp;/g;
      $line =~ s/</&lt;/g;
      $line =~ s/>/&gt;/g;
      print $line;
    }
    
    print '</PRE>';             # end preformatted output
  } else {
    print strong("<font color=red>Sorry, an error has occurred in reading the file \"$file\".</font>");
  }
} # end sub print_file
################################################################################

################################################################################
# Subroutine to READ-IN-FILES as user input on the screen (by uploading)
################################################################################
sub convert_to_unix {
  my ($newfile, $oldfile) = @_;

  if(open(OLDFILE, "$oldfile")) {
    open (NEWFILE, ">$newfile");

    while (<OLDFILE>) {
      s/\r//g;
      print NEWFILE;
    }

    close OLDFILE;
    close NEWFILE;
  } else {
    print strong("<font color=red>Sorry, an error has occurred in reading the file \"$oldfile\".</font>");
  }
} # end sub convert_to_unix
###############################################################################
