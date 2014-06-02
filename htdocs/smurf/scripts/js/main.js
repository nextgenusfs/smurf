/*
  Advanced Email Check
  By Website Abstraction (http://www.wsabstract.com) and
  Java-Scripts.net (http://www.java-scripts.net)
  Over 200+ free scripts here!
*/

var testresults=false

function checkemail(){
  var str=document.getElementById('email_val').value
  var filter=/^.+@.+\..{2,3}$/		
 
  if (filter.test(str))
    testresults=true
  else {
    alert("Please input a valid email address!")
    testresults=false
  }

  return (testresults)
}


var test=false;

function validator(registration){
  var email=document.getElementById('email_reg').value;
  var affiliation=document.getElementById('affiliation').value;
  var filter=/^.+@.+\..{2,3}$/;
	 
  if(email=="" && affiliation==""){
    alert("Please input a valid email address!\n\nPlease input an affiliation OR organization!");
    test=false;
  }
  else if(email != "" && affiliation == ""){
    alert("Please input a valid affiliation OR organization!");
    test=false;
  }
  else if(affiliation != "" && email == ""){
    alert("Please input a valid email address!");
    test=false;
  }
  else{
    test=true; 
    
    if(filter.test(email)){
      test=true;
    }
    else{
      alert("Please input a valid email address!");
      test=false;
    }

  }
  
  return(test);
}

function unhide(divID) {
  var id = document.getElementById(divID);
  id.style.display = (id.style.display == 'none' ? 'block' : 'none');
}

function upload_validator(upload){
  var FASTAfile = upload.FASTA.value;
  var COORDINATEfile = upload.coords.value;
  var FASTAbox = upload.FASTAseq.value;
  var COORDINATEbox = upload.coordinates.value;              
  var servergenomes = upload.genomes.options[upload.genomes.selectedIndex].value;
        
  if (FASTAfile != "" && COORDINATEfile != "" && FASTAbox != ""){
    alert("TOO MANY ARGUMENTS!\n\n\nUpload a protein FASTA or paste protein FASTA!\n\nUpload a coordinates file or paste coordinates!\n\nOR\n\nChoose genome input file from drop down menu below");
    test=false;
  } 
  else if (FASTAfile != "" && COORDINATEfile != "" && FASTAbox != "" && COORDINATEbox != ""){
    alert("TOO MANY ARGUMENTS!\n\n\nUpload a protein FASTA or paste protein FASTA!\n\nUpload a coordinates file or paste coordinates!\n\nOR\n\nChoose genome input file from drop down menu below");
    test=false;
  }
  else if (FASTAfile != "" && COORDINATEfile != "" && COORDINATEbox != ""){
    alert("TOO MANY ARGUMENTS!\n\n\nUpload a protein FASTA or paste protein FASTA!\n\nUpload a coordinates file or paste coordinates!\n\nOR\n\nChoose genome input file from drop down menu below");
    test=false;
  }
  else if (FASTAbox != "" && COORDINATEfile != "" && COORDINATEbox != ""){
    alert("TOO MANY ARGUMENTS!\n\n\nUpload a protein FASTA or paste protein FASTA!\n\nUpload a coordinates file or paste coordinates!\n\nOR\n\nChoose genome input file from drop down menu below");
    test=false;
  }
  else if ((servergenomes != "" && FASTAfile != "") || (servergenomes != "" && FASTAbox != "") || (servergenomes != "" && COORDINATEfile != "") || (servergenomes != "" && COORDINATEbox != "")){
    alert("TOO MANY ARGUMENTS!\n\n\nUpload a protein FASTA or paste protein FASTA!\n\nUpload a coordinates file or paste coordinates!\n\nOR\n\nChoose genome input file from drop down menu below");
    test=false;
  }
  else if(FASTAbox != "" && COORDINATEfile == "" && COORDINATEbox == ""){
    alert("Upload a coordinates file!\n\nOR\n\nPaste coordinates!");
    test=false;
  }
  else if(COORDINATEbox != "" && FASTAfile == "" && FASTAbox == ""){
    alert("Upload a protein FASTA!\n\nOR\n\nPaste protein FASTA!");
    test=false;
  }
  else if(FASTAbox != "" && COORDINATEbox != ""){
    test=true;
  }
  else if (FASTAfile=="" && COORDINATEfile =="" && FASTAbox == "" && COORDINATEbox == "" && servergenomes == ""){
    alert("Upload a protein FASTA or paste protein FASTA!\n\nUpload a coordinates file or paste coordinates!\n\nOR\n\nChoose genome input file from drop down menu below");
    test=false;
  }

  else if(FASTAfile != "" && COORDINATEfile=="" && COORDINATEbox==""){
    alert("Upload a coordinates file!\n\nOR\n\nPaste coordinates!");
    test=false;
  } 
  else if(COORDINATEfile != "" && FASTAfile=="" && FASTAbox==""){
    alert("Upload a protein FASTA!\n\nOR\n\nPaste protein FASTA!");
    test=false;
  }   
  else{
    test=true;
  }

  return(test);
}