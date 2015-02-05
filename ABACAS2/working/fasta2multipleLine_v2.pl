#!/usr/bin/perl

# Martin Hunt 11/03/09
# mh12@sanger.ac.uk

# takes a one-line-per-sequence fasta and converts to the evil multi-line format

use strict;
use warnings;

if ( @ARGV != 3){
    print "usage: fasta2multipleLine.pl <input fasta file> <output fasta file> <line width>\n";
	print "fasta has to be single lined\n" ;
    exit;
}

my $infile     = $ARGV[0];
my $outfile    = $ARGV[1];
my $line_width = $ARGV[2];

my $out_string = "";
my $tmp_string;

my %contig_seq = () ;
my @contig_names = () ;


print "reading contigs\n" ;
open (IN, "$infile") or die "oops!\n" ;
while (<IN>) {

    # print "$_" ;

    if (/^>(\S+)/) {

	my $seq_name = $1 ;
	my $seq = <IN> ;
	chomp($seq) ;

	$contig_seq{$seq_name} = $seq ;
	push(@contig_names, $seq_name) ;
    }
	
}
close(IN) ;


open F, ">", $outfile or die "Error opening $outfile";


foreach my $contig_name (@contig_names) {

	my $seq = $contig_seq{$contig_name} ;

	print "doing $contig_name\n" ;

	print F ">$contig_name\n" ;

	my $seq_len = length($seq) ;
	my $times = int($seq_len / $line_width) ;
	my $remainder = $seq_len % $line_width ;

	print "$seq_len $times\n" ;

	my $coords = 0 ;
	for (my $i = 0 ; $i < $times ; $i++ ) {
            my $out_string .= substr($seq,$coords,$line_width) . "\n";
		$coords += $line_width ;
            print F "$out_string" ;
        }
	if ($remainder != 0 ) {
        	my $out_string .= substr($seq,$coords,$line_width) . "\n";
        	print F "$out_string" ;
	}

}


close F;



