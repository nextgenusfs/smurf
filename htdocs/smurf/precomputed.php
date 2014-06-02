<?php
  require_once('/usr/local/common/web'.$_SERVER['WEBTIER'].'/templates/class.template.php');
  require_once('scripts/php/config.php');
  
  $template = new template();

  $side_menu[(!isset($_COOKIE['smurfs']) ? 2 : 3)]['class'] = 'no_sub side_active';
  $page_header = 'Precomputed Results';
  
  $content = '
    <p>For most seuquenced fungal species, predicted "backbone" genes encoding nonribosomal peptide synthases (NRPSs), polyketide synthases (PKSs), NRPS-like enzymes, PKS-like enzymes, Hybrid enzymes and demethylallyl tryptophan synthase (DMAT) and corresponding clusters can be downloaded <a href="ftp://ftp.jcvi.org/pub/software/smurf/">here.</a> </p> 
    <p><a href="javascript:unhide(\'genomes\');">Results are available for these sequenced fungal genomes:</a></p>
    <div id="genomes" style="display: none;">
      <table class="contenttable">
        <tbody>
          <tr class="tableHeader">       
            <td rowspan="1"><p>Sequenced Fungal Genomes</p></td>
          </tr>
          <tr class="tableRowEven">
            <td>
              <p>Aspergillus clavatus</p>
              <p>Aspergillus flavus</p>
              <p>Aspergillus fumigatus Af293</p>
              <p>Aspergillus fumigatus A1163</p>
              <p>Aspergillus nidulans</p>
              <p>Aspergillus niger CBS 513.88</p>
              <p>Aspergillus oryzae</p>
              <p>Aspergillus terreus</p>
              <p>Chaetomium globosum</p>
              <p>Coccidioides posadasii</p>
              <p>Coccidioides immitis</p>
              <p>Cryptococcus neoformans</p>
              <p>Fusarium graminearum</p>
              <p>Fusarium oxysporum</p>
              <p>Fusarium verticillioides</p>
              <p>Magnaporthe grisea</p>
              <p>Nectria haematococca</p>
              <p>Neosartorya fischeri</p>
              <p>Neurospora crassa</p>
              <p>Penicillium marneffei</p>
              <p>Phanerochaete chrysosporium</p>
              <p>Scerotinia sclerotiorum</p>
              <p>Stagonospora nodorum</p>
              <p>Trichoderma reesei</p>
              <p>Talaromyces stipitatus</p>
              <p>Uncinocarpus reesii</p>
              <p>Ustilago maydis</p>
            </td>
          </tr>
        </tbody>
      </table>
    </div>';
    
  $breadcrumb = array(
    array(
      'link' => 'precomputed.php',
      'name' => 'Precomputed',
    ),
  );

 
  $template->assign('javascript', array('scripts/js/main.js'));
  $template->assign('stylesheets', array('stylesheets/main.css'));
  $template->assign('title', $title.'Precomputed');
  $template->assign('top_menu', $top_menu);
  $template->assign('side_menu', $side_menu);
  $template->assign('page_header', $page_header);
  $template->assign('home_page', $home_page);  
  $template->assign('project_name', $project_name);
  $template->assign('main_content', $content);
  $template->assign('right_content', $right_content);
  $template->assign('site', 'SMURF');
  $template->assign('section_header', $section_header);
  $template->assign('breadcrumb', $breadcrumb);
  
  $template->display('3_column_fixed_width.tpl');
?>
