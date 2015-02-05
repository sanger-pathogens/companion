use strict;


### amout for overlap in Abacas2
my $overlap=100;


my %h;
open (F, shift) or die "Couldn't open file: $! ";
my $fasta=shift;
while (<F>){
	my @ar=split(/\t/);
	# 603052  606253  7       3200    3202    3194    99.44   946663  9599    1       1       DD2.chr2        NODE_3731_length_9535_cov_12.829471     
	
	my $chr=$ar[11];
	my $end=($ar[1]-$overlap);if ($end <1){$end=1};
	if (($ar[1]-$ar[0])>400) {
	  
	  foreach (($ar[0]+$overlap)..$end){
		$h{$chr}[($_-1)]=1;	
	  }
	}
	
  }

close(F);

my $ref_Seq=getFasta($fasta);

my $res;
my $found=0;
foreach my $chr (keys %h){
#	print "work on $chr\n";
	for my $pos (0..(length($$ref_Seq{$chr})-1)){
		if ($found && !defined($h{$chr}[$pos])){
			
			$res.=substr($$ref_Seq{$chr},$pos,1)
		}	
		elsif( !$found && !defined($h{$chr}[$pos])){
#			print "Start new stuff...\n";
			$res.="\n>fill.$chr.$pos\n".substr($$ref_Seq{$chr},$pos,1);
			$found=1;
		}
		else {
		   $found=0;	
		}
	}	
}

print $res;

sub getFasta{
	my $name = shift;
	
	open F, $name or die "no fasta to load $name: $!\n";
	my %h;
	my $fasta;
	while (<F>){
		chomp;
		if (/^>(\S+)/){
			$fasta=$1;	
		}
		else {
			$h{$fasta}.=$_
		}
	}	
	close(F);

	return \%h
}
