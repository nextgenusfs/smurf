<?php
  require_once('/usr/local/common/web'.$_SERVER['WEBTIER'].'/templates/class.template.php');
  require_once('scripts/php/config.php');
  
  $template = new template();

  $side_menu[(!isset($_COOKIE['smurfs']) ? 4 : 5)]['class'] = 'no_sub side_active';
  
  $content = '
    <p class="header2">Sequencing Centers:</p>
    <p><a href="http://www.jcvi.org/">J. Craig Venter Institute (JCVI)</a> (former TIGR)</p>
    <p><a href="http://www.broad.mit.edu/annotation/fgi/">The Broad Institute</a></p>			
    <p><a href="http://genome.jgi-psf.org/">Joint Genome Institute</a> (JGI)</p>			
    <p><a href="http://www.aist.go.jp/">AIST</a></p>
									
    <p class="header2">Fungal Genome Warehouses:</p>
    <p><a href="http://www.ncbi.nlm.nih.gov/projects/WGS/WGSprojectlist.cgi">GenBank WGS Projects</a></p>
    <p><a href="http://mips.gsf.de/projects/fungi">Munich Information Center for Protein Sequences</a> (MIPS)</p>
    <p><a href="http://www.cadre-genomes.org.uk/">Central Aspergillus Data Repository</a> (CADRE)</p>
    <p><a href="http://www.e-fungi.org.uk">E-fungi</a></p> 		
						
    <p class="header2">Other Fungal Resources:</p>
    <p><a href="http://www.aspergillus.org.uk/">The Aspergillus Website</a></p>
    <p><a href="http://www.aspergillusflavus.org/">AspergillusFlavus.org</a></p>
    <p><a href="http://fungalgenomes.org/">Fungal Genome Research</a></p>
    <p><a href="http://fungal.genome.duke.edu/">Resources for Fungal Comparative Genomics</a></p>
    <p><a href="http://gene.genetics.uga.edu/">Fungal Genome Resource</a> (at the University of Georgia)</p>
			
    <p class="header2">Protein Classification Databases:</p>
    <p><a href="http://pfam.sanger.ac.uk/">PFAM</a></p>
    <p><a href="http://www.tigr.org/TIGRFAMs/">TIGRFAM</a></p>
    <p><a href="http://www.ncbi.nlm.nih.gov/COG/">Clusters of Orthologous Groups</a> (COG/KOG)</p>
    <p><a href="http://www.tigr.org/jravel/nrps/">Prokaryotic PKS/NRPS analysis web site</a></p>
			
    <p class="header2">Pathway Databases:</p>
    <p><a href="http://www.genome.ad.jp/kegg/pathway.html">KEGG</a></p>
    <p><a href="http://www.reactome.org/">Reactome</a></p>
    <p><a href="http://metacyc.org/">MetaCyc</a></p>
    <p><a href="http://www.genego.com/">GeneGo</a></p>';  
    
  $breadcrumb = array(
    array(
      'link' => 'links.php',
      'name' => 'Links',
    ),
  );
  
  $template->assign('title', $title.'Links');
  $template->assign('top_menu', $top_menu);
  $template->assign('side_menu', $side_menu);
  $template->assign('page_header', 'Links');
  $template->assign('home_page', $home_page);  
  $template->assign('project_name', $project_name);
  $template->assign('main_content', $content);
  $template->assign('right_content', $right_content);
  $template->assign('site', 'SMURF');
  $template->assign('section_header', $section_header);
  $template->assign('breadcrumb', $breadcrumb);
  
  $template->display('3_column_fixed_width.tpl');
?>
