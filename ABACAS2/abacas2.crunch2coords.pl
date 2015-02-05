#! /usr/bin/perl -w
#
# File: assemlby.little.firsthit.pl
# Time-stamp: <20-Jul-2011 16:09:00 tdo>
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
  if (!defined($h{$ar[4]}) && $ar[0] > 40 ) {
	$h{$ar[4]}=1;
#	82 35 76 435 NODE_589_length_434_cov_17.797235 1048680 1049039
#	MAL10 none 
	
	my $strand=1;
	
	if ($ar[2] > $ar[3] && $ar[5] > $ar[6]) {
	  my $tmp=$ar[3];  $ar[3]= $ar[2];$ar[2]=$tmp;
	  $tmp=$ar[6];  $ar[6]= $ar[5];$ar[5]=$tmp;
	  
	  
	  $strand=1;
	}
	elsif ( $ar[5] > $ar[6]) {
	  $strand=-1;
	  my $tmp=$ar[3];  $ar[3]= $ar[2];$ar[2]=$tmp;
	  $tmp=$ar[6];  $ar[6]= $ar[5];$ar[5]=$tmp;
	}
	elsif ($ar[2] > $ar[3]) {
	  $strand=-1;
	}
	print "$ar[5]\t$ar[6]\t$ar[2]\t$ar[3]\t".($ar[6]-$ar[5])."\t".(abs$ar[3]-$ar[2])."\t$ar[1]\t1000\t1000\t1\t$strand\t$ar[7]\t$ar[4]\n"
  }
  
}
