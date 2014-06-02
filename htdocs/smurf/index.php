<?php
  require_once('/usr/local/common/web'.$_SERVER['WEBTIER'].'/templates/class.template.php');
  require_once('scripts/php/config.php');
  
  $template = new template();

  $side_menu[0]['class'] = 'no_sub side_active';
  
  $content = '
    <p class="header1">About SMURF</p>
    <p>Secondary Metabolite Unique Regions Finder is a web-based tool that finds secondary metabolite biosynthesis genes and pathways in fungal genomes. The predictions are based on <a href="http://pfam.sanger.ac.uk/">PFAM</a> and <a href="http://www.tigr.org/TIGRFAMs/">TIGRFAM</a> domain content as well as on a gene\'s chromosomal position. Precomputed clusters for most sequenced fungal genomes are available <a href="precomputed.php">here</a>. The software is described in Nora Khaldi, Fayaz T. Seifuddin, Geoff Turner, Daniel Haft, Ken Wolfe, William C. Nierman, Natalie D. Fedorova, 2008.</p>';
 
  $template->assign('title', $title.'About');
  $template->assign('top_menu', $top_menu);
  $template->assign('side_menu', $side_menu);
  $template->assign('page_header', $page_header);
  $template->assign('home_page', $home_page);  
  $template->assign('project_name', $project_name);
  $template->assign('main_content', $content);
  $template->assign('right_content', $right_content);
  $template->assign('site', 'SMURF');
  $template->assign('section_header', $section_header);
  
  $template->display('3_column_fixed_width.tpl');
?>
