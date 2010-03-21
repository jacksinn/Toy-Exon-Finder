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

#open contig file to read genetic sequence into variable $S
$handle = fopen($contigFile, 'r') or die("SUCK");
while(!feof($handle)){
    $line = trim(fgets($handle));
    #trim($line);
    $S = $S . $line;
}
fclose($handle);


?>