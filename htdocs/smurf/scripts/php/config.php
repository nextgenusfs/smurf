<?php
  $side_menu = array(
    array (
      'class' => 'no_sub',
      'link' => 'index.php',
      'menu_name' => 'About SMURF',
    ),
    array (
      'class' => 'no_sub',
      'link' => 'run_smurf.php',
      'menu_name' => 'Run SMURF',
    ),
  );
  
  // add upload menu page if logged in
  if(isset($_COOKIE['smurfs'])) { 
    $side_menu[] = array(
      'class' => 'subA',
      'link' => 'upload.php',
      'menu_name' => 'Upload Sequences',
    );  
  } // end if(isset($_COOKIE['smurfs']))

  $end_menu = array(
    array (
      'class' => 'no_sub',
      'link' => 'precomputed.php',
      'menu_name' => 'Precomputed',
    ),
    array (
      'class' => 'no_sub',
      'link' => 'faq.php',
      'menu_name' => 'FAQ',
    ),
    array (
      'class' => 'no_sub',
      'link' => 'links.php',
      'menu_name' => 'Links',
    ),
  );
  
  $side_menu = array_merge($side_menu, $end_menu);
  
  $right_content = array(
    array (
      'header' => 'NEW AND RETURNING USERS',
      'content' => '<p>Click <a href="run_smurf.php">here</a> to get started. </p>',
    ),
  );

  $page_header = 'Smurf - Secondary Metabolite Unique Regions Finder';
  $title = 'SMURF: ';
  $home_page = '/smurf';

  $section_header = array(
    'alttext' => 'Section Banner',
    'path' => 'images/banner-smurf.jpg',
  );
?>