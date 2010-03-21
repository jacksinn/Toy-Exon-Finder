<?php
########################################################################
#	Author: Steven C Jackson
#
#	Title:	CS 8550 Homework 2
#
#	Description: Toy Exon Finder
#
#	Use:	pass command line parameter of filename as first argument
#			e.g.
#			perl hw2.pl contig.txt
#			or /home/steven/./hw2.pl /home/steven/Documents/contig.pl
########################################################################

#variable setup
#echo $argv[1];
$contigFile = $argv[1];
$S = "";
$t = .68;
#echo $contigFile;
#open contig file to read genetic sequence into variable $S
$fh = fopen($contigFile, 'r') or die("SUCK");
while(!feof($fh)){
  $line = preg_replace('/ /', '', trim(fgets($fh)));
  #trim($line);
  #preg_replace("/ /", "//g", $line);
  $S = $S . $line;
}
#echo $S . "\n";
fclose($fh);

#print_r(getAllSignals($S));
#print_r(getAllOrfs($S));
#call TOYSCAN WHICH CALLS ALL THE OTHER FUNCTIONS
#RANDORF($S);
TOYSCAN_1($S, $t);


########################################################################
#FUNCTIONS
########################################################################
#getAllSignals takes a string of genetic code and parses out all items
#that may be signals and returns them in a hashed array
function getAllSignals($S){
	#$S = $_[0];
	$L = strlen($S);
	$psi = array();
	for($i=0; $i <= $L-2; $i++){
		$s = substr($S, $i, 2);
		if($s == "GT" || $s == "AG"){ array_push($psi, array($i => $s)); }
		if($i < $L-2){
			$s = substr($S, $i, 3);
			if($s == "TAG" || $s == "TGA" || $s == "TAA"){ $s = "TAG"; }
			if($s == "TAG" || $s == "ATG"){
				array_push($psi, array($i => $s));
			}
		}
	}
	return $psi;
}

function getAllOrfs($S){
	$psi = getAllSignals($S);
	$n = count($psi);
  $omega = array();
  
	for($i = 0; $i <= $n-2; $i++){
		$x1 = key($psi[$i]);  #position
    $s1 = $psi[$i][$x1];  #value
    
		if($s1 == "ATG" || $s1 == "AG"){
			if($s1 == "AG"){ $x1 = $x1 + 2; }
			for($j=$i+1; $j <= $n-1; $j++){
        $x2 = key($psi[$j]); #position;
        $s2 = $psi[$j][$x2]; #value
        if($s2 == "GT" || $s2 == "TAG"){
          if($s2 == "TAG"){
            $x2 = $x2 + 3;
            if($s1 == "ATG" && ($x2-$x1)%3 != 0){
              next;
            }
          }
          if($x1 < $x2){
            $rho = array($s1, $x1, $s2, $x2-1);
            if(hasOpenReadingFrame($rho, $S)){
              array_push($omega, $rho);
            }
          }
        }
			}
		}
	}
	return $omega;
	#echo "COUNT: " . $count . "\n";
}

function hasOpenReadingFrame($rho, $S){
  $s1 = $rho[0];
  $b  = $rho[1];
  $s2 = $rho[2];
  $e  = $rho[3];

  $delta = array(0, 1, 2);

	if($s1 == "ATG"){ $delta = array(0); }
	if($s2 == "TAG"){
		$e = $e - 3;
		$delta = array( ($e - $b + 1) %3 );
	}

	foreach($delta as $d){
		$stop = false;

		for($x = $b+$d; $x <= $e-2; $x = $x+3){
			$codon = substr($S, $x, 3);
			if($codon == "TAG" || $codon == "TGA" || $codon == "TAA"){
				$stop = true;
			}
		}
		if(!$stop){ return true; }
	}
	return false;
}

#Determining whether two ORFs overlap
function orfsOverlap($rho1, $rho2){
	$sb1 = $rho1[0];
	$b1 = $rho1[1];
	$se1 = $rho1[2];
	$e1 = $rho1[3];

	$sb2 = $rho2[0];
	$b2 = $rho2[1];
	$se2 = $rho2[2];
	$e2 = $rho2[3];

	return($b1 <= $e2 && $b2 <= $e1);
}

function RANDORF($S){
	$omega = getAllOrfs($S);
  #echo "OLD OMEGA::\n";
  #print_r($omega);
  $rho1POS = 0; #for unsetting array
  $rho2POS = 0; #for unsetting array
	foreach($omega as $rho1){
   // echo "O1::\n";
    //print_r($o1);
		foreach($omega as $rho2){
			if(orfsOverlap($rho1, $rho2)){
				if(rand(0, 100) < 50){
          //echo "LESS THAN 50\n";
          #echo $o1 . "\n";
					#pop $o1 off of @ohm
          unset($omega[$rho1POS]);
				}else{
          //echo "MORE THAN 50\n";
          #echo $o2 . "\n";
					#pop $o2 off of @ohm
          unset($omega[$rho2POS]);
				}
			}
      $rho2POS++;
		}
    $rho1POS++;
	}
  #echo "NEW OMEGA\n";
  #print_r($omega);
	return $omega;
}

function orfGcDensity($rho, $S){
	$s1 = $rho[0];
	$b  = $rho[1];
	$s2 = $rho[2];
	$e  = $rho[3];

	$gc = 0;
	$acgt = 0;

	for($i = $b; $i <= $e; $i++){
		$char = substr($S, $i, 1);
		if($char == "G" || $char == "C"){	$gc = $gc + 1;	}
		if($char == "A" || $char == "C" || $char == "G" || $char == "T"){
			$acgt = $acgt + 1;
		}
	}
	$gcdensity = $gc / $acgt;
  return $gcdensity;
}

function TOYSCAN_1($S, $t){
	$omega = getAllOrfs($S);
  $rhoCount = 0;
	foreach ($omega as $rho){
		if (orfGcDensity($rho, $S) < $t){
			unset($omega[$rhoCount]);
		}
    $rhoCount++;
	}

  $rho1POS = 0; #for unsetting array
  $rho2POS = 0; #for unsetting array
	foreach($omega as $rho1){
   // echo "O1::\n";
    //print_r($o1);
		foreach($omega as $rho2){
			if(orfsOverlap($rho1, $rho2)){
				if(orfGcDensity($rho1, $S) < orfGcDensity($rho2, $S)){
          //echo "LESS THAN 50\n";
          #echo $o1 . "\n";
					#pop $o1 off of @ohm
          unset($omega[$rho1POS]);
				}else{
          //echo "MORE THAN 50\n";
          #echo $o2 . "\n";
					#pop $o2 off of @ohm
          unset($omega[$rho2POS]);
				}
			}
      $rho2POS++;
		}
    $rho1POS++;
	}
  echo "NEW OMEGA\n";
  print_r($omega);
	return $omega;
}
?>