#! /usr/bin/perl -w
#
# File: aba2.sort.coords.pl
# Time-stamp: <18-Nov-2009 17:56:03 tdo>
# $Id: $
#
# Copyright (C) 2009 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
#

use strict;
use Data::Dumper;


my $file=shift;
my $contigFile = shift;
my $step = shift;

#nucmer 'n' or promer 'p'
my $method=shift;

if (!defined($step)) {
  $step=0;
}

my $minAlignmentLength=500;


if (!defined($method) || $method eq 'nucmer') {
  $method='n';
}
elsif ($method eq 'promer') {
  print "Using Promer ... \n";
  $minAlignmentLength=50;
} 
else {
  die "Please specify method nucmer / promer\n";
  
  
}

if (!defined($contigFile)) {
  die "coords fiel> <contigFile> <step 0 - just order 1 - split>\n";
  
}



open (F,$file) or die "die gimme file.";
my @F=<F>;
close(F);

my %h;
my %contigLength;
my $prePos;
my $preCo='';
my $preStr='';
my $preRef='';
my $firstLine;

my %h2;
my %ContigsTokeep;

my ($refS,$refE,$contigS,$contigE,$alignRef,$alignContig,$id,$similarity,$refLength,$contigLength,$strand1,$strand2,$reference,$contig,$dummy,$dummy2);

# load all alignments in hash
foreach (@F) {
  if (/^\d/) {

	if ($method eq 'n') {
	  ($refS,$refE,$contigS,$contigE,$alignRef,$alignContig,$id,$refLength,$contigLength,$strand1,$strand2,$reference,$contig)
	  = split(/\s+/);
	}
	else {
	  ($refS,$refE,$contigS,$contigE,$alignRef,$alignContig,$id,$similarity,$dummy,$refLength,$contigLength,$strand1,$strand2,$reference,$contig)
	  = split(/\s+/);
	}
	
	$h{$contig}{$alignContig}=$_;
	$contigLength{$contig}=$contigLength;
	if ($strand2==1) {
	  # positive
	  $h2{$contig}{$contigS}="$reference,$alignContig,$strand2,$contigS,$contigE,$refS";
	}
	else {
	  #contig is reersed, start und stop also inversted
	  $h2{$contig}{$contigE}="$reference,$alignContig,$strand2,$contigE,$contigS,$refS";
	}
	
  }
}

my @ar;
my $countCut=0;


# write alginment by order of alignmentoverlap to array @ar;
foreach my $contig (sort keys %h) {
  foreach my $align (sort {$b <=> $a} keys %{$h{$contig}}) {
	push @ar, $h{$contig}{$align};
  }
}




my $res;
my %FirstLine;

my %contigs2Analyse;

foreach (@ar) {
  if (/^\d/) {
  
	if ($method eq 'n') {
	  ($refS,$refE,$contigS,$contigE,$alignRef,$alignContig,$id,$refLength,$contigLength,$strand1,$strand2,$reference,$contig)
	  = split(/\s+/);
	}
	else {
	  ($refS,$refE,$contigS,$contigE,$alignRef,$alignContig,$id,$similarity,$dummy,$refLength,$contigLength,$strand1,$strand2,$reference,$contig)
	  = split(/\s+/);
	}
	
	# save the best it
	if ($preCo ne $contig) {
	  $FirstLine{$contig}=$_;
	  
	  $firstLine=$_;
	  $preRef=$reference;
	  $prePos=$refS;
	  $preCo=$contig;
	  $preStr=$strand2;	  
	  $res.=$_;
	  $ContigsTokeep{$contig}=1;
	}
	else {
	  if ($preRef ne $reference ||
#		  $preStr != $strand2   ||
		  (abs($prePos-$refS)>(2*$contigLength))
		 ) {
		$countCut++;
		
		if ($alignContig>$minAlignmentLength) {
		  $contigs2Analyse{$contig}=1;
		  
		}
	  }
	  else {
		$res.=$_;
	  }
	}
	
  }
}

if ($step == 0 ) {
  print "Save the new coords file\n";
  
  open (F,"> $file.ordered") or die "problem to write $file.ordered \n";
  print F $res;
  close(F);
  print "$countCut contigs alignemt thrown out\n";

  print "Contig to splice: ", scalar keys  %contigs2Analyse ,"\n";
  foreach my $contig (keys %contigs2Analyse) {
	print "$contig\n";
	
  }

  exit;
  
}



my $count=0;
my $preContigMappingEnd;
my %cutContig;

foreach my $contig ( keys  %contigs2Analyse ){
  my $cutStart=0;
  my $contigLength=$contigLength{$contig};
  my $countit=0;
  $preRef='';
  $prePos='';
  $preStr=0;	  
  $preContigMappingEnd=-1;
  #  print Dumper  $h2{$contig};
  
  
  $countit=0;
  
  foreach my $pos (sort {$a <=> $b} keys %{ $h2{$contig} } ){
	
	my ($ref,$alignLength,$strand,$contigStart,$contigEnd,$refS)=split(/,/,$h2{$contig}{$pos});
	
#	print ">>> Position ($contig): $pos - $countit : $contigStart - $contigEnd strand: $strand alignLength $alignLength\ncut start: $cutStart, preContigMapping: $preContigMappingEnd\n$ref,$alignLength,$strand,$contigStart,$contigEnd,$refS\n\n";
	
	# first initialization
	if (! $countit) {
	  $preRef=$ref;
	  $prePos=$refS;
	  $preStr=$strand;	  
	  $preContigMappingEnd=$contigEnd;
#	  print "Set $contigEnd $preContigMappingEnd\n";
	  
	}
	elsif ($preRef ne $ref ||
		   $preStr != $strand   ||
		   (abs($prePos-$refS)>(2*$contigLength))
		  ) {
	  if ($alignLength>500) {
		
		push @{$cutContig{$contig}}, "$cutStart,$preContigMappingEnd";
#		print "put $contig ($cutStart,$contigStart,$preContigMappingEnd) - $strand \n";
		
		$preRef=$ref;
		$prePos=$refS;
		$preStr=$strand;	  
		$preContigMappingEnd=$contigEnd;
		$cutStart=$contigStart;
		
	  }
	  else {
#		print "ignoe parte $contig ($cutStart,$contigStart,$preContigMappingEnd) $strand \n";
	  }
	}
	
	
	# same contig strand, update information
	else {
	  $preRef=$ref;
	  $prePos=$refS;
	  $preStr=$strand;
	  
	  $preContigMappingEnd=$contigEnd;
	}
	$countit++
  }
  
  push @{$cutContig{$contig}}, "$cutStart,".($contigLength-1);
#  print Dumper $cutContig{$contig};
#  system ("grep $contig All.morphed.coords")
	
	
}
 
my $ref_fasta=loadFasta($contigFile);

$res='';


foreach my $contig (keys %cutContig) {
  my $num=1;
  foreach my $bit (@{$cutContig{$contig}}) {
	print "Get bit $bit as Part number $num\n";
	my ($from,$to)=split(/,/,$bit);
	$res.=">$contig.Part$num $from:$to ".($to-$from+1) ."\n";
	$res.=substr($$ref_fasta{$contig},$from,($to-$from+1))."\n";
	$num++;
  }
}

foreach my $contig (sort keys %ContigsTokeep) {
  if (!defined($cutContig{$contig})) {
  	$res.=">$contig\n$$ref_fasta{$contig}\n";
  }
}


open(F,"> $contigFile.splitted") or die "probelms $!\n";
print F $res;

close(F);

print "$countCut contigs alignemt thrown out\n";
print "$count contigs are going to be doubled\n";






sub loadFasta{
  my $fasta=shift;
  
  open F, $fasta or die "prob couldn't find fasta file: $fasta $!  \n";

  my %h;
  my @ar;
  my $name;
  
  while (<F>) {
        chomp;
        if (/^>(\S+)/){
		  $name=$1;
		}
		else {
		  chomp;
		  $h{$name}.=$_;
		}
	  }
  return \%h;
}
