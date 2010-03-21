<?php
   /*
    #!/usr/bin/perl
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

use strict;
use warnings;

sub getAllOrfs;
sub getAllSignals;
sub hasOpenReadingFrame;
sub RANDORF;
sub orfGcDensity;

#variable setup
my $S = "";
my $t = .68;

#open contig file to read genetic sequence into variable $S
open FILE, $ARGV[0] or die $!;
while(my $line = <FILE>){
	chomp($line);
	$S = $S . $line;
}
close FILE;


#getAllOrfs($S);
TOYSCAN_1($S, $t);


########################################################################
#FUNCTIONS
########################################################################
sub getAllOrfs{
	my $S = $_[0];
	my @omega = ();
	my %psi = getAllSignals($S);
	my @sorted_psi = sort keys %psi; #to maintain same order
	#size of hash
	my $n = keys(%psi);
	#while(($index, $signal) = each(@sorted_psi)){
	#	print $index . " : " . $signal . "\n";
	#}
	#foreach $key (@sorted_psi){
	#	print $key .  " | " . $value . "\n";
	#}
	my $count = 0;
	for(my $i = 0; $i <= $n-2; $i++){
		my $x1 = $sorted_psi[$i]; #position
		my $s1 = $psi{$sorted_psi[$i]}; #value
		#($s1, $x1) = $sorted_psi[$i];
		print $s1 . " | " . $x1 . "\n";

		if($s1 eq "ATG" or $s1 eq "AG"){
			if($s1 eq "AG"){
				#print "AG whatwhat\n";
				$x1 = $x1 + 2;
			}
			for(my $j=$i+1; $j <= $n-1; $j++){
					my $x2 = $sorted_psi[$j]; #position;
					my $s2 = $psi{$sorted_psi[$j]}; #value
					#print $s2 . "************\n";
					if($s2 eq "GT" or $s2 eq "TAG"){
						if($s2 eq "TAG"){
							$x2 = $x2 + 3;
							if($s1 eq "ATG" and ($x2-$x1)%3 != 0){
								next;
							}
						}
						if($x1 < $x2){
							my @rho = ($s1, $x1, $s2, $x2-1);
							#print "^^^^^^^" . $s2 . "WOOT\n";
							$count++;
							if(hasOpenReadingFrame(@rho, $S) == 1){
								#print "======YAY======\n";
								push(@omega, @rho);
							}
						}
					}
			}
		}
	}
	return @omega;
	#print "COUNT: " . $count . "\n";
}

#getAllSignals takes a string of genetic code and parses out all items
#that may be signals and returns them in a hashed array
sub getAllSignals{
	my $S = $_[0];
	my $L = length($S);
	my %psi;
	for(my $i=0; $i <= $L-2; $i++){
		my $s = substr($S, $i, 2);

		if($s eq "GT" or $s eq "AG"){
			#PUSH ONTO LIST (PSI)
			#associative arrays use curly braces >.<
			$psi{$i} = $s;
		}
		if($i < $L-2){
			$s = substr($S, $i, 3);
			if($s eq "TAG" or $s eq "TGA" or $s eq "TAA"){
		    	$s = "TAG"
			}
			if($s eq "TAG" or $s eq "ATG"){
				#PUSH ONTO LIST PSI
				$psi{$i} = $s;
			}
		}
	}
	return %psi;
}


sub hasOpenReadingFrame{
	my @rho = ($_[0], $_[1], $_[2], $_[3]);
	my $S = $_[4];

	#foreach(@rho){print "\t$_\n---\n";}

	my $s1 = $rho[0];
	my $b = $rho[1];
	my $s2 = $rho[2];
	my $e = $rho[3];

	#print "-------------" . $s2 . "\n";
	my @delta = (0, 1, 2);

	if($s1 eq "ATG"){
		@delta = (0);
	}
	if($s2 eq "TAG"){
		$e = $e - 3;
		@delta = ( ($e - $b + 1) %3 );
	}
	foreach my $d (@delta){
		my $stop = 0;

		for(my $x = $b+$d; $x <= $e-2; $x = $x+3){

			my $codon = substr($S, $x, 3);
			#print "!!! " . $x . " : " . $x+2 . ": $codon\n";
			if($codon eq "TAG" or $codon eq "TGA" or $codon eq "TAA"){
				$stop = 1;
			}
		}
		if($stop == 0){ return 1; }
	}
	return 0;
}

#Determining whether two ORFs overlap
sub orfsOverlap{
	my @rho1 = $_[0];
	my @rho2 = $_[1];

	my $sb1 = $rho1[0];
	my $b1 = $rho1[1];
	my $se1 = $rho1[2];
	my $e1 = $rho1[3];

	my $sb2 = $rho2[0];
	my $b2 = $rho2[1];
	my $se2 = $rho2[2];
	my $e2 = $rho2[3];

	if($b1 <= $e2 and $b2 <= $e1){
		return 1;
	}else{
		return 0;
	}
}

sub RANDORF{
	my $S = $_[0];

	my @ohm = getAllOrfs($S);

	foreach my $o1 (@ohm){
		#print "O1: $o1\n";
		foreach my $o2 (@ohm){
			if(orfsOverlap($o1, $o2) == 1){
				if(int(rand(100)) < 50){
					#pop $o1 off of @ohm
				}else{
					#pop $o2 off of @ohm
				}
			}
		}
	}

	return @ohm;
}

sub orfGcDensity{
	my @rho = ($_[0], $_[1], $_[2], $_[3]);
	my $S = $_[4];

	my $s1 = $rho[0];
	my $b = $rho[1];
	my $s2 = $rho[2];
	my $e = $rho[3];

	my $gc = 0;
	my $acgt = 0;

	for(my $i = $b; $i <= $e; $i++){
		my $char = substr($S, $i, 1);
		if($char eq "G" or $char eq "C"){
			$gc = $gc +1;
		}
		if($char eq "A" or $char eq "C" or $char eq "G" or $char eq "T"){
			$acgt = $acgt + 1;
		}
	}
	my $gcdensity = $gc / $acgt;
}

sub TOYSCAN_1{
	my $S = $_[0];
	my $t = $_[1];

	my @ohm = getAllOrfs($S);

	foreach my $rho (@ohm){
		if (orfGcDensity($rho, $S) < $t){
			#pop rho off of ohm
		}
	}
	foreach my $rho1 (@ohm){
		foreach my $rho2 (@ohm){
			if(orfsOverlap($rho1, $rho2) == 1){
				if(orfGcDensity($rho1, $S) < orfGcDensity($rho2, $S)){
					#pop rho1 off of ohm
				}else{
					#pop rho2 off of ohm
				}
			}
		}
	}
	return @ohm;
}
    */
?>