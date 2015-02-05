#! /usr/bin/perl -w
#
# File: abacas2.gene.reference.pl
# Time-stamp: <13-Jan-2011 20:59:37 tdo>
# $Id: $
#
# Copyright (C) 2011 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description: This program will generate a reference just from the
# CDS.  
#

use strict;


my $embl=shift;
my $resultName=shift;

my $seq=embl2fasta($embl);

my @new;

open F, "$embl " or die "Problem open $embl in Embl2Fasta:$! \n";

while (<F>) {
  chomp;
  if (/source/) {
  }
  elsif (/^FT   \S+\s+(\S+)$/) {
	my @ar=split(/,/);
	foreach (@ar) {
	  if ( /(\d+)\.\.(\d+)/){
		for my $i ($1..$2) {
		  $new[$i]=substr($seq,($i-1),1);
		}
	  }
	  
	}
  }
}

my $res;
for $_ (1..(scalar(@new)-1)) {
  if (defined($new[$_])) {
	$res.=$new[$_]
  }
  else {
	$res.="N"
  }
}
  my ($name)=$embl =~ /^(\S+)\.embl/;
print ">$name\n";

print $res;


sub embl2fasta{
  my $embl = shift;
  
  open F, "$embl " or die "Problem open $embl in Embl2Fasta:$! \n";
  my ($name)=$embl =~ /^(\S+)\.embl/;
  
  my $fasta;
  
  while (<F>) {
	if (/^SQ/) {
	  while (<F>) {
		if (/\/\//) {
		  
		  last;
		}
		## get away space and number of seq
		s/\d+//g;
		s/\s+//g;
		$fasta.=$_;
		
	  }
	}
  }
  
  return $fasta
}
