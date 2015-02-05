#!/usr/bin/perl -w
use lib '/home/otto/bin';
use strict;


### will cut a fasta in individual files with the name >(\S*) A$
# ARGV[1] will say how much sequences...
SeperateSequences($ARGV[0],$ARGV[1]);
sub SeperateSequences{
  my $FileName = shift;
  my $Amount = shift;

  if (!defined($Amount)) {
	$Amount = 99999999999;
  }
  open (FILE, $FileName) or die "Couldn't open file to Seqparate $FileName $!\n";
  
  my $Counter = 0;
  my $Out;
  my @NameFiles;
  while (<FILE>) {
	if ($Counter < $Amount) {
	  
	  if (/^>(\S+)/) {
		print $_,"\n";
		
		close(FILEOUT);
		my $Name = $1;#."_$Counter.fasta";
		open (FILEOUT, "> $Name") or die "Problem do $Name file  $!\n";
		push @NameFiles, $Name;
		$Counter++;
		print FILEOUT $_;
	  }
	  else {
		print FILEOUT $_;
	  }
	}
  }
  print "$Counter Sequences where generated\nDone.\n";
  return \@NameFiles;
}
