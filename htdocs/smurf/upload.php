<?php
  // redirect to run_smurf if not logged in
  if(!isset($_COOKIE['smurfs'])) { 
    header('Location: run_smurf.php');
  } // end if(!isset($_COOKIE['smurfs']))

  require_once('/usr/local/common/web'.$_SERVER['WEBTIER'].'/templates/class.template.php');
  require_once('scripts/php/config.php');
  
  $template = new template();

  $side_menu[2]['class'] = 'subA side_active';
  
  $content = '
    <p class="header1">Two types of files are required</p>
    <ul>
      <li class="no_sub">MultiFASTA file containing protein sequences. Headers must contain proteinID as the first element. <a href="javascript:unhide(\'MultiFASTA\');">View MultiFASTA file example</a></li>
    </ul>
    <div id="MultiFASTA" style="display: none">
      <table class="contenttable">
        <tbody>
          <tr class="tableHeader">       
            <td rowspan="1"><p>MultiFASTA file example</p></td>
          </tr>
          <tr class="tableRowEven">
            <td>
              <p>>AFUA_6G05350</p>
              <p>MSVRQVPKQRRTLIEFSAFKEVPCMLFCVAMFFGYIGFFNPIFYIEAFAIQKHAMGETLAFHL</p>
              <p>ISILNATSVPGRIVPGILGLRFGPLNILLGSAIISGILSLCWIAIYNAGPLIVLAVLYGFSGAFVS</p>
              <p>LLAVALTTLNLNLQTLRTRMGMCSLLCGFGSLCRAPVAGAVLDNTRSYLGVQLYSGLTIGTTGV</p>
              <p>LLFFANHLKRRTN*</p>
              <p>>AFUA_6G05320</p>
              <p>MHSWIDIPLRDGKDADEWELTESKAIANIDQWYQEGRLLFPRDSLQEVRDQLRKPLKKGDSI</p>
              <p>YVKGFDGSMYEWPVKTGNIKTHWEADNNTGEEKQVDWMLLVLEKSDIKVLEDPEGFLEMDDASLV</p>
              <p>FSCSPNICVKEIVVDIARPLALIWCTVKDKDPEHAL*</p> 
              <p>>AFUA_6G05310</p>
              <p>MTAIHNKSQKQSWNLSWNLSWNLSWNLSWNQSWNQSWNLSWNLSWNQSRKHNQEYKTMSGH</p>
              <p>RPAQKRKERNLLAHAKIYVFATIYLVDALREQCLKSLHRDLSNFVLNRQTITNVLDLLEY</p>
              <p>TYDHTGRQEPGGRCSLRTLVIHYISFENTRFRRILDDHGGMASDLVGEAC*</p> 
            </td>
          </tr>
        </tbody>
      </table>
    </div>
    <ul>
      <li class="no_sub">Tab delimited TXT file containing the following elements in the order shown below: <a href="javascript:unhide(\'genecoordinates\');">View gene coordinates file example</a></li>
    </ul>
    <ol>
      <li class="no_sub">protein ID</li>
      <li>chromosome/contig</li>
      <li>5\' gene start</li>
      <li>3\' gene stop</li>
      <li>protein name/function/definition (if available)</li>
    </ol>
    <div id="genecoordinates" style="display: none">
      <table class="contenttable">
        <tbody>
          <tr class="tableHeader">       
            <td rowspan="1"><p>protein ID</p></td>
            <td rowspan="1"><p>chromosome/contig</p></td>
            <td rowspan="1"><p>5\' gene start</p></td>
            <td rowspan="1"><p>3\' gene stop</p></td>
            <td rowspan="1"><p>protein name/function/definition (if available)</p></td>
          </tr>
          <tr class="tableRowEven">
            <td><p>AFUA_6G05350</p></td>
            <td><p>101</p></td>
            <td><p>13578</p></td>
            <td><p>11957</p></td>
            <td><p>aspartic-type endopeptidase (OpsB), putative</p></td>	
          </tr>
          <tr class="tableRowOdd">
            <td><p>AFUA_6G05320</p></td>
            <td><p>101</p></td>
            <td><p>18874</p></td>
            <td><p>20201</p></td>
            <td><p>purine nucleoside phosphorylase I, inosine and guanosine-specific</p></td>
          </tr>
          <tr class="tableRowEven">
            <td><p>AFUA_6G05310</p></td>
            <td><p>101</p></td>
            <td><p>21397</p></td>
            <td><p>20403</p></td>
            <td><p>nucleolus protein required for cell viability, putative</p></td>
          </tr>	
        </tbody>
      </table>
    </div>
    <form id="upload" onsubmit="return upload_validator(this)" action="/cgi-bin/smurf/upload.cgi" method="post" enctype="multipart/form-data">
      <p><b>Upload a protein FASTA file</b> <input type="file" name="FASTA" size="20" /></p>
      <p><b>Or, Paste protein FASTA</b> <textarea name="FASTAseq" cols="60" rows="5"></textarea></p>
      <p><b>Upload a coordinates file</b> <input type="file" name="coords" size="20" /></p>
      <p><b>Or, Paste coordinates</b> <textarea name="coordinates" cols="60" rows="5"></textarea></p>
      <p><input type="submit" value="Submit" /> <input type="reset" value="Reset" /></p>

      <p><b>If you would like to use processed input files instead, select one of the following genomes that are available on the server</b></p>
      <p><select name = "genomes">
        <option value=""> Choose a genome...</option>
        <option value="aclavatus"> A. clavatus</option>
        <option value="aflavus"> A. flavus</option>
        <option value="afumigatus"> A. fumigatus Af293</option>
        <option value="afumigatus_A1163"> A. fumigatus A1163</option>
        <option value="anidulans"> A. nidulans</option>
        <option value="aniger"> A. niger 513.88</option>
        <option value="aoryzae"> A. oryzae</option>
        <option value="aterreus"> A. terreus</option>
        <option value="cglobosum"> C. globosum</option>
        <option value="cimmitis"> C. immitis</option>
        <option value="cneoformans"> C. neoformans</option>
        <option value="cposadasii"> C. posadasii</option>
        <option value="fgraminearum"> F. graminearum</option>
        <option value="foxysporum"> F. oxysporum</option>
        <option value="fverticillioides"> F. verticillioides</option>
        <option value="mgrisea"> M. grisea</option>
        <option value="ncrassa"> N. crassa</option>
        <option value="nfischeri"> N. fischeri</option>
        <option value="pchrysogenum_newest"> P. chrysogenum</option>
        <option value="pchrysosporium"> P. chrysosporium</option>
        <option value="pmarneffei"> P. marneffei</option>
        <option value="sclerotiorum"> S. clerotiorum</option>
        <option value="snodorum"> S. nodorum</option>
        <option value="treesei"> T. reesei</option>
        <option value="tstipitatus"> T. stipitatus</option>
        <option value="umaydis"> U. maydis</option>
        <option value="ureesii"> U. reesii</option>
      </select></p>
      <p><input type="submit" value="Submit" onclick="window.location=document.upload.genomes.options[document.upload.genomes.selectedIndex].value" /> <input type="reset" value="Reset" /></p>
    </form>	
';
  
  $breadcrumb = array(
    array(
      'link' => 'run_smurf.php',
      'name' => 'Run SMURF',
    ),
    array(
      'link' => 'upload.php',
      'name' => 'Upload Sequences',
    ),
  );
  
  $template->assign('javascript', array('scripts/js/main.js'));
  $template->assign('title', $title.'Upload Sequences');
  $template->assign('top_menu', $top_menu);
  $template->assign('side_menu', $side_menu);
  $template->assign('page_header', 'Upload Sequences');
  $template->assign('home_page', $home_page);  
  $template->assign('project_name', $project_name);
  $template->assign('main_content', $content);
  $template->assign('right_content', $right_content);
  $template->assign('site', 'SMURF');
  $template->assign('section_header', $section_header);
  $template->assign('breadcrumb', $breadcrumb);
  
  $template->display('3_column_fixed_width.tpl');
?>
