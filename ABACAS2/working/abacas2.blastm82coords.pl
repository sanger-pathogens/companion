#! /usr/bin/perl -w
#
# File: assemlby.little.firsthit.pl
# Time-stamp: <29-Jul-2011 14:42:46 tdo>
# $Id: $
#
# Copyright (C) 2011 by Pathogene Group, Sanger Center
#
# Author: Thomas Otto
#
# Description:
#


use strict;


my %h;


while (<STDIN>) {
  my @ar=split(/\s+/);
  if (!defined($h{$ar[0]}) && $ar[2] > 40 ) {
	$h{$ar[4]}=1;
#	82 35 76 435 NODE_589_length_434_cov_17.797235 1048680 1049039
#	MAL10 none 
	
	my $strand=1;
	
	if ($ar[6] > $ar[7] && $ar[8] > $ar[9]) {
	  my $tmp=$ar[7];  $ar[7]= $ar[6];$ar[6]=$tmp;
	  $tmp=$ar[9];  $ar[9]= $ar[8];$ar[8]=$tmp;
	  
	  
	  $strand=1;
	}
	elsif ( $ar[8] > $ar[9]) {
	  $strand=-1;
	  my $tmp=$ar[9];  $ar[9]= $ar[8];$ar[8]=$tmp;
		
	}
	elsif ($ar[2] > $ar[3]) {
	  $strand=-1;
	}
	my $lgth1=1000;
	my $lgth2=1000;
	
	if (defined($ar[13])) {
	  $lgth2=$ar[12];
	  $lgth1=$ar[13];
	}
	
	print "$ar[8]\t$ar[9]\t$ar[6]\t$ar[7]\t".($ar[9]-$ar[8])."\t".(abs$ar[7]-$ar[6])."\t$ar[2]\t$lgth1\t$lgth2\t1\t$strand\t$ar[1]\t$ar[0]\n"
  }
  
}
