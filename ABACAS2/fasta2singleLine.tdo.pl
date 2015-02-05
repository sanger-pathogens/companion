use strict;

my %h;

my $count=0;
while (<STDIN>){
	if (/^>/){ 
		if ($count){
			print "\n"; 	
		}
		print ;
		$count++;
		}
	else {
		chomp;
		print ;	
	}	
	
}
