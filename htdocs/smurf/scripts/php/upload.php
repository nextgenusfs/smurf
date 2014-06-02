<?php
  echo '<pre>';
  debug('rm -rf /opt/www/smurf/tmp/59Nd');
  debug('mkdir /opt/www/smurf/tmp/59Nd');
  debug('chmod 777 /opt/www/smurf/tmp/59Nd');

  debug('hmmsearch --cut_tc /opt/www/smurf/data/hmmsPKSNRPS/AMP-binding.hmm /opt/www/smurf/data/outfiles/OUT_aclavatus.input > /opt/www/smurf/tmp/59Nd/hmm_AMP.out');
  debug('hmmsearch --cut_tc /opt/www/smurf/data/hmmsPKSNRPS/Acyl_transf_1.hmm  /opt/www/smurf/data/outfiles/OUT_aclavatus.input > /opt/www/smurf/tmp/59Nd/hmm_Acyl_transf.out');
  debug('hmmsearch --cut_tc /opt/www/smurf/data/hmmsPKSNRPS/Condensation.hmm /opt/www/smurf/data/outfiles/OUT_aclavatus.input > /opt/www/smurf/tmp/59Nd/hmm_Condensation.out');
  
  echo '</pre>';
  
  function debug($str) { 
    echo $str.' ('.system($str, $ret).' - '.$ret.')'."\n";
  } // end function debug($str)
?>