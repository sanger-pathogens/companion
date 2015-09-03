package TilingGraph;

# tdo 01.11.09: include that when new contig, the contig with max
# overhang left, is take



#To DO lkist
# inlcude %lefthandPosContig

use 5.008008;
use strict;
use warnings;
use Data::Dumper;


require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration       use Mapping::MappingTools ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw(

) ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

our @EXPORT = qw( startTiling

);

our $VERSION = '0.01';
my $DEBUG=1000;
my $MIN_OVERLAP=0.1;

my $CHECK_OVERLAP=$ENV{"ABA_CHECK_OVERLAP"};

# double chec this, clean
my $countIT=1000;
my $NAME;
my $ref_leftHandPosContig;

my $MAX_CONTIG_END=500;
my $MAX_GAP_SIZE=10000;
my $MIN_GAP_SIZE=100;
my $GAP_SIZE=100;


# Below is stub documentation for your module. You'd better edit it!

=head1 NAME

Protocols::Assembly::TillingRaph - Perl extension for performing
assemblies guided with a reference.

=head1 SYNOPSIS

use Protocols::Assembly::TillingGraph;

=head1 DESCRIPTION

This model will be able to join contigs of different assemblers to a
pseudo molecule. The ordering is done by mapping the congitig with
nucmer (mummer package) against a given Reference: doNumcer()
For better performance the contigs itself will be blasted against each
others to define there potential overlap, see
blastOverlaps.pl blast contigs.fasta contigs.fasta output_prefix
-report_overlaps -minLength 50 -endTolerance 10 -annotate
of mh12.

This porgram will load the nucmer and the overlap files, and define if
they overlap, in case they will be joined. if not, a given number of
nnn's will be included.
Other programs were tested:
show-tiling of the mummer package is buggy, and good mapping contigs
are ignored.
minimus2 didn't performe well and cannot use the order of the contig
by a reference.

=head1 EXPORT

startTiling

=head1 SEE ALSO

http://scratchy.internal.sanger.ac.uk/wiki/index.php/Team_133

# Preloaded methods go here.

=head1 FUNCTIONS

=head2 startTilling - This function will call all programs

        Usage           :

        Arg [1]         : LSF queue for jobs

        Arg [2]         : file of mapped lanes to be summarised

        Arg [3]         : faidx index file for reference

        Example         :

        Description     : Summarises the coverage of the mapped lanes specified in the lanes.fofn file

        Returntype      : none

        Author          : Jacqueline McQuillan  E<lt>jm15@sanger.ac.uk<gt>

=cut

sub startTiling
{
  my $showTilingFile   = shift;
  my $contigFile       = shift;
  my $resultname       = shift;


  if (!defined( $CHECK_OVERLAP)) {
	$CHECK_OVERLAP=1;
  }

  $countIT=0;
  $NAME=$resultname;


  ### load contigs
  my $ref_contigSeq= loadFasta($contigFile);


  ### load the overlp of the contigs, done see before
  my ($ref_overlapGraph,$ref_repeatContigs);  # =
                                              # loadOverlap($overlapFile);


  ### load the positionGraph against the reference
  my ($ref_positionGraph,$ref_contigLength,$ref_contigStrand) = loadShowTiling($showTilingFile);
  debug(10,"Print PostionGraph");


  ### construct the sequence
  my ($ref_Sequence,$ref_referenceCoverage,$ref_contigGFF,$ref_gapsGFF) = printPositionGraph($ref_positionGraph,$ref_contigLength,$ref_overlapGraph,$ref_contigSeq,$ref_contigStrand);


  ### WRITE THE FILES
  save($resultname,$ref_Sequence,$ref_referenceCoverage,$ref_contigGFF,$ref_gapsGFF);

}

sub save{
  my $resultName    = shift;
  my $ref_Sequence  = shift;
  my $ref_referenceCoverage = shift;
  my $ref_contigGFF = shift;
  my $ref_gapsGFF      = shift;
  my $count=0;


  foreach my $chr ( sort keys %$ref_Sequence) {


	open (F, "> $resultName.$chr.fna") or die "Couldn't open $resultName: $!\n";
	print F ">$chr\n";
	print F join('',$$ref_Sequence{$chr});
	print F "\n";
	close(F);

#	open (F, "> $resultName$chr.gaps.gff") or die "Couldn't open $resultName.gff: $!\n";
#	print F $$ref_gapsGFF;
#	close(F);

	open (F, "> $resultName.$chr.refCoverage.plot") or die "Couldn't open $resultName.$chr.refCoverage.plot: $!\n";
	### check for zeros:
	foreach (0..(scalar(@{$$ref_referenceCoverage{$chr}})-1)) {
	  if (!defined($$ref_referenceCoverage{$chr}[$_])) {
		$$ref_referenceCoverage{$chr}[$_]="0\n";
	  }
	}
	my $res = join('',@{$$ref_referenceCoverage{$chr}});
	print F $res;
	close(F);

	open (F, "> $resultName.$chr.contigs.gff") or die "Couldn't open $resultName.gff: $!\n";
	print F $$ref_contigGFF{$chr};
	close(F);
  }


}

=head2 makeOverlap

        Usage           :

        Arg [1]         : Overlap report

        Example         :

        Description     : generate the overlap of two contigs and returns  %contigOverlaps

        Returntype      : reference to %contigOverlaps

        Author          : Jacqueline McQuillan E<lt>jm15@sanger.ac.uk<gt>

=cut

sub makeOverlap{
  my $contig1 = shift;
  my $contig2 = shift;
  my $ref_contigSeq = shift;
  my $NAME= shift;


  $countIT++;
  system("rm tmp.$NAME.$countIT.* &> /dev/null");

  open (F, "> tmp.$NAME.$countIT.fasta" ) or die "Couldnprobelm s\n";
  print F ">$contig1\n".join('',@{$$ref_contigSeq{$contig1}}),"\n>$contig2\n", join('',@{$$ref_contigSeq{$contig2}}),"\n";
  close(F) or die "Problems, I couldn't close the file handle: $! \n";

  while (!(-R "tmp.$NAME.$countIT.fasta" && -s "tmp.$NAME.$countIT.fasta")) {
	sleep(3);

  }

 	sleep(3);
  ## extract the contigs
  my $call="export PERL5LIB; blastOverlaps.pl blast tmp.$NAME.$countIT.fasta tmp.$NAME.$countIT.fasta tmp.$NAME.$countIT.run -report_overlaps -minLength 50 -endTolerance $MAX_CONTIG_END -annotate -annotate_only";
 # print $call,"\n";
#  $ENV{PERL5LIB}="/software/badger/lib/perl5:/software/pathogen/psu_cvs/genlib/perl/src:/software/pathogen/external/lib/site_perl:/software/pathogen/external/lib/perl:/software/pathogen/external/lib/perl/lib:/software/pathogen/external/lib/perl/lib/site_perl:/software/pathogen/external/lib/perl5:/software/pubseq/PerlModules/Modules:/software/badger/lib/perl5:/software/noarch/badger/lib/perl5:/nfs/users/nfs_t/tdo/bin/oldPerl:/nfs/users/nfs_t/tdo/pathogen/user/mh12/perl/lib:/software/pathogen/projects/protocols/lib/perl5";
 # system("set ");
#  print("$call");


  !system("$call " ) or die "Problem with: $call";

  while (!(-R "tmp.$NAME.$countIT.run.blast.annotated.gz" && -s "tmp.$NAME.$countIT.run.blast.annotated.gz")) {
	sleep(3);
  }
  sleep(3);

  # system ("zcat  tmp.$NAME.$countIT.run.blast.annotated.gz ");
  my ($ref_overlapGraph,$ref_repeatContigs) = loadOverlap("tmp.$NAME.$countIT.run.blast.annotated.gz");

  #  system("rm tmp.$NAME.$countIT.*");

  return ($ref_overlapGraph,$ref_repeatContigs)
}
=head2 loadOverlap

        Usage           :

        Arg [1]         : Overlap report

        Example         :

        Description     : load the Overlap into a hash %contigOverlaps

        Returntype      : reference to %contigOverlaps

        Author          : Jacqueline McQuillan E<lt>jm15@sanger.ac.uk<gt>

=cut

sub loadOverlap
{
# quick if rep OVERLAP All.overlap.blast.annotated >
# All.overlap.blast.annotated.Overlap
  my $overlapFile = shift;


  my %contigOverlaps;
  my %repeatContigs;


  open(F, "gunzip -c $overlapFile | ") ||
	die "Couldn't open File $overlapFile: $!\n";

  # get file in array, just a bit quicker
  my @F=<F>;
  close(F);

  my ($contig1,$contig2,$id,$overlap,$start1,$end1,$start2,$end2);

  foreach (@F) {
	#	print ;

	if (/OVERLAP/) {
	  ($contig1,$contig2,$id,$overlap,$start1,$end1,$start2,$end2) = $_ =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+\d+\s+\d+\s(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s/;
	  #by definition, there can be just one overlap between contigs;
	  # but they are summetric!
	  #	  print "$id 90 \n";

	  if (($id>=100 && $overlap >= 20)
		  ||
		  (($id>=99 && $overlap >= 100))
		 ) {

		#construct the overlap adjactace graph. Each node is
		#(orientation contig1, ori. contig2,

		#get orientation: general code not just for blast
		my $contigOrientation1="+";
		my $contigOrientation2="+";
		if ($end1<$start1) {
		  $contigOrientation1="-";
		}
		if ($end2<$start2) {
		  $contigOrientation2="-";
		}
		$repeatContigs{$contig1}++;
		#	print "$contigOrientation1,$contigOrientation2,$start1,$end1,$start2,$end2";
		$contigOverlaps{$contig1}{$contig2}="$contigOrientation1,$contigOrientation2,$start1,$end1,$start2,$end2";
	  }
	}

  }

  if (!defined($contigOverlaps{$contig1}{$contig2})) {
	my %longestOverlap;
	$longestOverlap{$contig1}{$contig2}=0;
	$longestOverlap{$contig2}{$contig1}=0;

	foreach (@F) {
	  #	  print ;

	  ($contig1,$contig2,$id,$overlap,$start1,$end1,$start2,$end2) = $_ =~ /^(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+\d+\s+\d+\s(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s/;
	  #by definition, there can be just one overlap between contigs;
	  # but they are summetric!
	  #	  print "$id 90 \n";


	  if ($id>=99 && $overlap >= 500 && $contig1 ne $contig2)		{

#		print "took it!!!\n";

		#construct the overlap adjactace graph. Each node is
		#(orientation contig1, ori. contig2,

		#get orientation: general code not just for blast
		my $contigOrientation1="+";
		my $contigOrientation2="+";
		if ($end1<$start1) {
		  $contigOrientation1="-";
		}
		if ($end2<$start2) {
		  $contigOrientation2="-";
		}
		if (!defined($longestOverlap{$contig1}{$contig2}) || $overlap > $longestOverlap{$contig1}{$contig2}) {
		  #		  print "Booked it";

		  #	print "$contigOrientation1,$contigOrientation2,$start1,$end1,$start2,$end2";
		  $contigOverlaps{$contig1}{$contig2}="$contigOrientation1,$contigOrientation2,$start1,$end1,$start2,$end2";
		  $longestOverlap{$contig1}{$contig2}=$overlap
		}

	  }
	}
  }
  #  print Dumper %contigOverlaps;

  return (\%contigOverlaps,\%repeatContigs)
}

=head2 loadShowTiling

        Usage           :

        Arg [1]         : Showtiling arquive of all cointigs agianst the reference. IMPORTANT to generate it with with "show-coords -bcloTq -I 95 ALL.Chab.filter.delta > ALL.Chab.filter.coords". The delta-filta should have be done with delta-filter -q. There will be a functino doing this.

        Description     : Load the order of the contigs in a hash(chromsom) of array postion hash which contigs covere

        Returntype      : reference of this structure.

        Author          : Jacqueline McQuillan E<lt>jm15@sanger.ac.uk<gt>

=cut

sub loadShowTiling
{
  my $showTilingFile = shift;

  open(F, $showTilingFile) ||
	die "Couldn't open File $showTilingFile: $!\n";

  # get file in array, just a bit quicker
  my @F=<F>;
  close(F);

  my %positionGraph;
  my %contigMappingLength;
  my %contigLength;
  my %contigStrand;

  my $ref_positionGraph;

  my @overlapPositions=('','','','','','');
  my $previousContig='';
  my $previousReference;
  my ($refS,$refE,$contigS,$contigE,$alignRef,$alignContig,$id,$similarity,$refLength,$contigLength,$strand1,$strand2,$reference,$contig,$dummper);
  my $countStrand=0;
  my $iteration=0;

  foreach (@F) {
	if (/^\d+/) {
	  $iteration++;

	  # parse through the show-coords 50641 58342 1 7794 7702 7794
	  #498452 34679 1.55 22.47 Morphed.Chab.10.chr02 Contig_0000760
	  #242546  245015  19586   17135   2470    2452    98.74   1168589 20540   1       -1      Morphed.Chab.10.chr07   contig2142

	  if (defined($ENV{'ABA_COMPARISON'}) && $ENV{'ABA_COMPARISON'} eq "promer") {
	  ($refS,$refE,$contigS,$contigE,$alignRef,$alignContig,$id,$similarity,$dummper,$refLength,$contigLength,$strand1,$strand2,$reference,$contig)	= split(/\s+/);

	  }
	  else {
		($refS,$refE,$contigS,$contigE,$alignRef,$alignContig,$id,$refLength,$contigLength,$strand1,$strand2,$reference,$contig)
		= split(/\s+/);
	  }


	  # longer the hits, more it counts. $strand is 1 / -1


	  if ($previousContig ne $contig && $iteration > 1 ) {

		if ($previousContig ne '' ) {

		  $ref_positionGraph=fillPositionGraph($ref_positionGraph,\@overlapPositions,$previousReference,$previousContig,\%contigLength,\%contigMappingLength);


		}

		$$ref_leftHandPosContig{$contig}=9999999;
		@overlapPositions=($refS,$refE,$contigS,$contigE,$strand1,$strand2);
		print "1>> $iteration >$contig - Countstrand $countStrand $strand2 $alignRef \n";
		if ($countStrand >= 0) {
		  $contigStrand{$previousContig}=1
		}
		else {
		  $contigStrand{$previousContig}=-1;
		}


		$contigLength{$contig}=$contigLength;
		$contigMappingLength{$contig}+=$alignContig;
		$previousContig=$contig;
		$previousReference=$reference;

		$countStrand=0

	  } # end $previousContig ne $contig

	  elsif ($iteration > 1){
		if ($previousReference eq $reference    && #$strand2 eq $contigStrand{$contig} &&
			abs($refS-$overlapPositions[0]) < 2*   $contigLength{$contig}
		   ) {

		  $contigMappingLength{$contig}+=$alignContig;
		  #check if Reference is mapping reverse
		  if ($refS<$overlapPositions[0]) {
			$overlapPositions[0]=$refS
		  }
		  if ($overlapPositions[1] < $refE) {
			$overlapPositions[1]=$refE
		  }


		  # check also for the contigs, which could be a problem with
		  # show-coords -r
		  if ($contigS<$overlapPositions[2]) {
			$overlapPositions[2]=$contigS
		  } else {
			$overlapPositions[3]=$contigE
		  }
		}
		else {
		  debug(40,"ignore line $_");

		}


	  } ## in first iteration
	  else {
		$previousContig=$contig;
		$previousReference=$reference;
		$contigMappingLength{$contig}+=$alignContig;
		$contigLength{$contig}=$contigLength;
		$contigMappingLength{$contig}+=$alignContig;
		@overlapPositions=($refS,$refE,$contigS,$contigE,$strand1,$strand2);

	  }

	  $countStrand+=($strand2*$alignRef);
	  print "2>> $iteration >$contig - Countstrand $countStrand $strand2 $alignRef \n";
	} ## end if /d


	# Fill last line;
  }

  ### last contigs...
  $ref_positionGraph=fillPositionGraph($ref_positionGraph,\@overlapPositions,$reference,$contig,\%contigLength,\%contigMappingLength);
  @overlapPositions=('','','','','',''); # jsut in case ofor empty line;

  print "1>> $iteration >$contig - Countstrand $countStrand $strand2 $alignRef \n";
  if ($countStrand >= 0) {
	$contigStrand{$contig}=1
  }
  else {
	$contigStrand{$contig}=-1;
  }

  $contigLength{$contig}=$contigLength;
  $contigMappingLength{$contig}+=$alignContig;

  return ($ref_positionGraph,\%contigLength,\%contigStrand)
}


=head2

        Usage           : Private

        Description     : Fills the Postion graph

        Returntype      : none

        Author          : Jacqueline McQuillan E<lt>jm15@sanger.ac.uk<gt>

=cut

sub fillPositionGraph
{
  my $ref_positionGraph      = shift;
  my $ref_positions          = shift;
  my $referenceName          = shift;
  my $contigName             = shift;
  my $ref_contigLength       = shift;
  my $ref_contigMappinglength= shift;

  # here check could be made over the contig legnth ofr the Contigh
  # Mapping Length!!!

  debug(120,"$referenceName $contigName ".$$ref_positions[0] ." - ". $$ref_positions[1]);

  my $overlap=($$ref_positions[1]-$$ref_positions[0]);

  if (defined($$ref_contigLength{$contigName}) && ($overlap/$$ref_contigLength{$contigName}*100>$MIN_OVERLAP)) {

	if (defined($$ref_leftHandPosContig{$contigName}) && $$ref_positions[0] > $$ref_leftHandPosContig{$contigName}) {
	  $$ref_leftHandPosContig{$contigName}= $$ref_positions[0]
	}

	debug(40,"Fill $contigName from $$ref_positions[0] to $$ref_positions[1]");

	for ($$ref_positions[0]..$$ref_positions[1])
	  {
		$$ref_positionGraph{$referenceName}[$_]{$contigName}=$contigName;
	  }

  }

  else {
	debug(120,"Contig $contigName not engouh  overlap $overlap / " .$$ref_contigLength{$contigName} );

  }

  return $ref_positionGraph
}
=head2

        Usage           :

        Arg [1]         : LSF queue for jobs

        Arg [2]         : file of mapped lanes to be summarised

        Arg [3]         : faidx index file for reference

        Example         :

        Description     : write the poostios of the contigs, in like abacas order

        Returntype      : none

        Author          : Jacqueline McQuillan E<lt>jm15@sanger.ac.uk<gt>

=cut

sub printPositionGraph
{
  my $ref_positionGraph  = shift;
  my $ref_contigLength   = shift;
  my $ref_overlapGraph   = shift;
  my $ref_contigSeq      = shift;
  my $ref_contigStrand   = shift;

  my %visitedContigs;
  my %layout;
  my %sequence;
  my $ref_sequence = \%sequence;

  my %referenceCoverage;

  my %core;
  my $setcore=0;



  if (-f "list.core.txt") {
	open (F,"list.core.txt") or die "Weird error openeing list.core.txt\n";
	while (<F>) {
	  if (/^(\S+):\s+(\d+)\.+(\d+)/) {
		$core{$1}[0]=$2;
		$core{$1}[1]=$3;
	  }
	}
	$setcore=1;

  }

#  print Dumper %$ref_positionGraph;

  ### the gff annotation files for the contigs and the gaps
  my ($ref_contigGFF,$ref_gapsGFF);


  # loop over chromosomes of reference
  foreach my $chr (keys %$ref_positionGraph) {
	debug(10,"Work on $chr");

	# set loop varialbes per chromosome
	my $amountUsedContigs=0;     # maount of build in contigs
	my $amountGaps       = 0;  # amount of gaps
	my $contigStartSeq   = 1;  # where the contigs starts to overlap to
                               # the next
	my $gapSizeCount     = 0;
	my $pseudoPosition   = 1;


	my $chrLength=(scalar(@{$$ref_positionGraph{$chr} })-1);

	my $pos       =1;    # actual position of chrosmome $chr
	my $actualContig='';
	my (@actualPos)=('','');



	#### MAIN LOOP
	while ($pos <= $chrLength) {

	  if ($DEBUG>80) {
		if (($pos%2000)==0) {

		}
	  }
	  ### get there the coverage of the reference
	  if (!defined($$ref_positionGraph{$chr}[$pos])){
		$referenceCoverage{$chr}[$pos].="0\n";
	  }
	  else {
		$referenceCoverage{$chr}[$pos].= keys %{$$ref_positionGraph{$chr}[$pos]};
		$referenceCoverage{$chr}[$pos].="\n";;
	  }

	  if ( !($setcore) or  # if core.list is given just use range
		   (
			defined($core{$chr}[0]) and # set a range for a chr
			$core{$chr}[0] <= $pos and
			$pos <= $core{$chr}[1]

		   )
		 ) {
	  	if (($pos%2000)==0) {
#		print "pos $pos\n";

		}

		### no contig at all on position $pos
		if (!defined($$ref_positionGraph{$chr}[$pos]) &&
			$actualContig eq '') {
		  ### nothing found no contig in sight.  count the gap, amount
		  ### of n's
		  $gapSizeCount++
		}
		### again a contig, after gap
		elsif ($actualContig eq '') {
		  $actualContig=getContigNewContig($$ref_positionGraph{$chr}[$pos],$ref_contigLength);
		  debug(40, "Found contig $actualContig on pos $pos STRAND  $$ref_contigStrand{$actualContig}");
		  (@actualPos)=($pos,$pos);
		}

		### actualcontig is coming to an end
		elsif (! defined($$ref_positionGraph{$chr}[$pos]{$actualContig})
			  ){
		  debug (40,"Search ($pos) overlap Contig $actualContig (".$$ref_contigLength{$actualContig}.") \t$actualPos[0] \t $actualPos[1] ");



		  my ($newContig,$posNew)=getNextContig($$ref_positionGraph{$chr},$ref_contigLength,$pos);
		  debug(100," Acutal $actualContig NewContig - $newContig STRAND  $$ref_contigStrand{$actualContig};");

		  ($ref_sequence,$actualContig,$contigStartSeq,$amountGaps,$pseudoPosition,$ref_contigGFF)
			=buildOverlap($ref_sequence,
						  $actualContig,$newContig,$ref_contigSeq,
						  $ref_contigStrand,
						  $ref_overlapGraph,$chr,
						  $contigStartSeq,$amountGaps,
						  $pos,
						  $pseudoPosition,$ref_contigGFF,
						  $gapSizeCount,$actualPos[0]);

		  # new gap might start, put count to zero
		  $gapSizeCount=0;

		  if ($actualContig ne '') {
			$actualPos[0]=$pos;  $actualPos[1]=$pos;
			$pos=($posNew-1); # the position might be new set from

			# getNextContig
#			$pos+=10;

		  }
		  else {
			debug(100,"gaps");

			# do this, as a contig is found, so this will be used.
		#	$pos+=100;

		  }

		}

		### just contig on the same contig
		else {
		  $actualPos[1]=$pos;
		}

	  } # end if core

	  #increment the position in thePositionGraph
	  $pos++;

	} # end while pos
	### END MAINLOOP

	if (defined($actualContig) && $actualContig ne '') {
	  ($ref_sequence,$actualContig,$contigStartSeq,$amountGaps,$pseudoPosition,$ref_contigGFF)
		  =buildOverlap($ref_sequence,
						$actualContig,"",$ref_contigSeq,
						$ref_contigStrand,
						$ref_overlapGraph,$chr,
						$contigStartSeq,$amountGaps,
						$pos,
						$pseudoPosition,$ref_contigGFF);
}
	$amountGaps--;
	debug(0,"Statistics:\nAmount of gaps $chr: $amountGaps\nUsed contigs: $amountUsedContigs\n");

  } # end loop chromoms


  return ($ref_sequence,\%referenceCoverage,$ref_contigGFF,$ref_gapsGFF);

}

sub getContig{
  my $ref_Contigs      = shift;
  my $ref_ContigLength = shift;

  my $length=0;
  my $longestContig='';
  foreach my $contig ( keys %$ref_Contigs ) {
	if ($$ref_ContigLength{$contig} > $length){
	  $length=$$ref_ContigLength{$contig};
	  $longestContig=$contig;
		}
  }
  return $longestContig;
}

sub getContigNewContig{
  my $ref_Contigs      = shift;
  my $ref_ContigLength = shift;

  my $minPos=99999999999;
  my $longestContig='';
  foreach my $contig ( keys %$ref_Contigs ) {
	if ($$ref_leftHandPosContig{$contig} < $minPos){
	  $minPos=$$ref_leftHandPosContig{$contig};
	  $longestContig=$contig;
	}
  }
  return $longestContig;
}


sub getContigBorder{
  my $ref_Contigs      = shift;
  my $ref_ContigLength = shift;
  my $pos              = shift;
  my $offset           = shift;

  #pos is already for the new cointg
  my $previousPos=($pos-$offset);


  my $length=0;
  my $longestContig='';
#  print "$pos -- $previousPos\n";

#  print Dumper $$ref_Contigs[$pos];
#  print Dumper $$ref_Contigs[$previousPos];

  foreach my $contig ( keys %{$$ref_Contigs[($pos)]} ) {
	# check if contig is left and right of gap and is the longs
	if (defined($$ref_Contigs[($pos-$offset)]{$contig}) && # position
                                                           # in
                                                           # $offset pos
		defined($$ref_Contigs[$previousPos]{$contig}) && # also before
                                                         # in the graph
		$$ref_ContigLength{$contig} > $length){
	  $length=$$ref_ContigLength{$contig};
	  $longestContig=$contig;
	}
  }
  return $longestContig;
}
sub getNextContig{
  my $ref_Contigs      = shift;
  my $ref_ContigLength = shift;
  my $pos              = shift;



  my $maxWalk=1000;

  my $notFound=1;

  # check, if short gap overlapping contig exists.
  my $longestContig=getContigBorder($ref_Contigs,$ref_ContigLength,$pos,500);
  print "longest $longestContig \n";

  if ($longestContig eq '') {
	$longestContig=getContigBorder($ref_Contigs,$ref_ContigLength,$pos,75);
	if ($longestContig eq '') {
	  $longestContig=getContigBorder($ref_Contigs,$ref_ContigLength,$pos,10);
	  if ($longestContig eq '') {
		$pos+=10;
		while ($pos < ($pos+$maxWalk) && $notFound) {
		  if (defined($$ref_Contigs[$pos])) {
			$notFound=0;

			$longestContig=getContig($$ref_Contigs[$pos],$ref_ContigLength);

		  }
		  $pos++;
		}
	  }
	}

  }
  return($longestContig,($pos-1))
}


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
		  my @a=split(//);
		  foreach (@a) {
			push @{ $h{$name}}, $_;
		  }
		}
	  }
  return \%h;
}

=head2

        Usage           : debug(level,msg)

        Arg [1]         : level when to printsub

        Arg [2]         : msg as string

       Example         : debug(100,"this ist an error")

        Description     : Prinst debug informat

        Returntype      : none

        Author          : Jacqueline McQuillan E<lt>jm15@sanger.ac.uk<gt>

=cut

sub debug
{

  my $code     = shift;
  my $msg      = shift;

  if ($code <= $DEBUG) {
	print $msg."\n";

  }
#code for the function goes here

}
=head2

        Usage           : privagter

        Arg [1]         : LSF queue for jobs

        Arg [2]         : file of mapped lanes to be summarised

        Arg [3]         : faidx index file for reference

        Example         :

        Description     : Will join two contigs

        Returntype      : string hash, Startposition new contigs

        Author          : Jacqueline McQuillan E<lt>jm15@sanger.ac.uk<gt>

=cut

sub buildOverlap
{
  my $ref_sequence     = shift;
  my $actualContig     = shift;
  my $newContig        = shift;
  my $ref_contigSeq    = shift;
  my $ref_contigStrand = shift;
  my $ref_overlapGraph_= shift;
  my $chr              = shift;
  my $contigStartSeq   = shift;
  my $amountGaps       = shift;
  my $position         = shift;
  my $pseudoPosition   = shift;
  my $ref_contigGFF    = shift;
  my $gapSizeCount     = shift;
  my $start_Ref_map    = shift;


  my $contigUsedLength=0;

#  print Dumper $$ref_contigSeq{$actualContig};

  my $contigLength=(scalar(@{$$ref_contigSeq{$actualContig}}));

  my $ref_cont=\@{$$ref_contigSeq{$actualContig}};

  # is coming from a gap
  if ($gapSizeCount>0) {
	debug(23, "Gap away: $chr - $position - $gapSizeCount ".($contigLength-$position)." -  $start_Ref_map ");
	my $gapEnd=($start_Ref_map-1);
	my $gapStart=( $gapEnd-$gapSizeCount);


$gapSizeCount=$GAP_SIZE;
	### tdo Nov 2014 - gap size fixed to 100bp
	### tdo Sep 2014 maximal gap size
	if ($gapSizeCount> $MAX_GAP_SIZE) {
	  $gapSizeCount=$GAP_SIZE;

	}
# taken out, just one gap, of $GAP_SIZE
#	$$ref_contigGFF{$chr}.="unknown\tGAP\tGAP\t$pseudoPosition\t".($pseudoPosition+$gapSizeCount-1)."\t0\t+\t.\tcolor=9;label=\"GAP\";note=\"Gap+of+size+$gapSizeCount+REFERENCE=$chr:$gapStart-$gapEnd\"\n";

#	$$ref_sequence{$chr} .= makeN($gapSizeCount);
#	$pseudoPosition+=$gapSizeCount;
  }

  my ($ref_overlapGraph,$ref_nothing);
  print "Building Overlap $actualContig $newContig\n";

  if ($CHECK_OVERLAP && defined($newContig) && $newContig ne '' ) {
	($ref_overlapGraph,$ref_nothing) = makeOverlap($actualContig,$newContig,$ref_contigSeq,$chr);
  }
#  print Dumper $ref_overlapGraph;

  ### here we have to also check, if the real end of the contig is
  ### doing the overlap.
  # (a) if contig Strand actual Contig +, end must overlap
  my $overlapOK=0;

  if ($CHECK_OVERLAP && defined($newContig) && defined($$ref_overlapGraph{$actualContig}{$newContig})) {
	$overlapOK=checkOverlap($ref_overlapGraph,$actualContig,$newContig,$ref_contigStrand,$ref_contigSeq)
  }

  if ($overlapOK) {
	debug(50,"Report overlap $actualContig $newContig fro $chr");
	# my ($ref_overlapData) =
	# split(/\,/,$$ref_overlapGraph{$actualContig}{$newContig});

	my (@overlapData) = split(/\,/,$$ref_overlapGraph{$actualContig}{$newContig});


	if ($DEBUG>500) {
	  system("grep $actualContig.*$newContig *.overlap");
	  system("egrep \"$actualContig|$newContig\" *.coords");
	}

	# join the actual contig to the pseudo sequence
	if ($$ref_contigStrand{$actualContig}==1){
	  debug(50,"Report +  overlap $actualContig $newContig $contigStartSeq-1 $overlapData[3] ");
	  $contigUsedLength=(($overlapData[2])-($contigStartSeq)+1);

	  $$ref_sequence{$chr} .= join('',@$ref_cont[($contigStartSeq-1)..($overlapData[2]-2)]);
	  $$ref_contigGFF{$chr}.="unknown\tContig\tContig\t$pseudoPosition\t".($pseudoPosition+$contigUsedLength)."\t0\t+\t.\tcolor=3;label=\"contig=$actualContig\";note=\"$actualContig+($contigLength)\"\n";;
	  if ($overlapData[1] eq"+") {
		$contigStartSeq=$overlapData[4]
	  }
	  else {
		$contigStartSeq=((scalar(@{$$ref_contigSeq{$newContig}})-1)-$overlapData[4]+1);
	  }
	}
	else {
	  # here the $contigSart is pretty high
	  $contigUsedLength=($contigLength-$overlapData[3]-($contigStartSeq-1));
	  debug (23," $overlapData[3] : $contigLength : $contigStartSeq \n");

	  $$ref_sequence{$chr} .=revcomp( join('',@$ref_cont[($overlapData[3])..($contigLength-$contigStartSeq)]));
	  $$ref_contigGFF{$chr}.="unknown\tContig\tContig\t$pseudoPosition\t".($pseudoPosition+$contigUsedLength-1)."\t0\t-\t.\tcolor=3;label=\"REVERSED.contig=$actualContig\";note=\"$actualContig+($contigLength)\"\n";
	  debug(50,"Report -  overlap $actualContig with $newContig \nTake sequencd from ".($overlapData[3])." to  $contigLength minus ".($contigStartSeq-1)." length $contigUsedLength. New StartSeq  $overlapData[5] ");
	  if ($overlapData[1] eq "-") {
		$contigStartSeq=$overlapData[5]
	  }
	  else {
		$contigStartSeq=((scalar(@{$$ref_contigSeq{$newContig}})-1)-$overlapData[5]);
	  }

	}

	# set pos where to set now feature
	$pseudoPosition+=($contigUsedLength-1);

	$$ref_sequence{$chr} .="\n";



	$actualContig=$newContig;

  }    ###

	### no overlap, just put the contig
  else {
	debug(40,"Position $position");

	$amountGaps++;
	$contigUsedLength=($contigLength-($contigStartSeq-1));
	debug(40,"Gap found: use $contigUsedLength of $contigLength. $actualContig ".($contigStartSeq)." ".(scalar( @$ref_cont) -1));


	if ($$ref_contigStrand{$actualContig}==1){
	  $$ref_sequence{$chr} .=join('',@$ref_cont[($contigStartSeq-1)..($contigLength-1)]);
	  $$ref_contigGFF{$chr}.="unknown\tContig\tContig\t$pseudoPosition\t".($pseudoPosition+$contigUsedLength-1)."\t0\t+\t.\tcolor=4;label=\"contig=$actualContig\";note=\"$actualContig+($contigLength)\"\n";;
	}
	else {
	  $$ref_sequence{$chr} .=revcomp(join('',@$ref_cont[0..($contigLength-1-($contigStartSeq-1))]));

	  $$ref_contigGFF{$chr}.="unknown\tContig\tContig\t$pseudoPosition\t".($pseudoPosition+$contigUsedLength-1)."\t0\t-\t.\tcolor=4;label=\"REVERSED.contig=$actualContig\";note=\"$actualContig+($contigLength)\"\n";
	}

	# adapt pos with contigs size
	$pseudoPosition+=($contigUsedLength);

	$$ref_contigGFF{$chr}.="unknown\tGAP\tGAP\t".($pseudoPosition)."\t".($pseudoPosition+$GAP_SIZE-1)."\t0\t+\t.\tcolor=7;label=\"GAP2\";note=\"Gap+as+no+Overlap+found.+Size+prefixed+$MIN_GAP_SIZE\"\n";

	$$ref_sequence{$chr} .= makeN($GAP_SIZE);
	$pseudoPosition+=($GAP_SIZE);
	print 	$pseudoPosition."\n";



	#include here the gpas GFFnd the contigs GFF.
	$contigStartSeq=1;
	$actualContig='';

  }

  return ($ref_sequence,$actualContig,$contigStartSeq,$amountGaps,$pseudoPosition,$ref_contigGFF)
}


#
#	$overlapOK=checkOverlap($ref_overlapGraph,$actualContig,$newContig,$ref_contigStrand,$ref_contigSeq)

sub checkOverlap{
  my ($ref_overlapGraph,$actualContig,$newContig,$ref_contigStrand,$ref_contigSeq) = @_;

  print "Check overlap on  $actualContig $$ref_contigStrand{$actualContig} and $newContig $$ref_contigStrand{$newContig}\n";
  my $overlapOK=1;
  my $lenActual=(scalar(@{$$ref_contigSeq{$actualContig}})-1);
  my $lenNew=(scalar(@{$$ref_contigSeq{$newContig}})-1);
  my $strandActual=$$ref_contigStrand{$actualContig};
  my $strandNew   =$$ref_contigStrand{$newContig};
  print "$lenActual $lenNew \n";

  #Position, where the contigs is doing the alignment
  my $alignmentActual='';
  my $alignmentNew='';
#          'contig00172' => {
#                             'NODE_134_length_517_cov_135.433273' => '+,+,1,48,500,547'
#                           },

  # where 1 is actual contig and 2 is newContig
  my ($strand1,$strand2,$start1,$end1,$start2,$end2)=
	split(/,/,$$ref_overlapGraph{$actualContig}{$newContig});

  if (
	  ($strand1 eq '+')
	 ){
	if ($start1 < $MAX_CONTIG_END)  {
	  $alignmentActual='Begin';
	  debug(10,"Actual Overlap + Begin")
	}
	elsif ($end1 > ($lenActual-$MAX_CONTIG_END)) {
	  $alignmentActual='End';
	  debug(10,"Actual Overlap + End")
	}
	else {
	  #	  $overlapOK=0;
	  debug(10,"Deeep shit Actual +")
	}
  }
  else {
	die "this shouldn't occur!\n";


  }

  print "Strand Actual $strandActual $alignmentActual\n";

  if (
	  ($strandActual == 1 && $alignmentActual eq 'Begin')
	  or
	  ($strandActual eq "-1" && $alignmentActual eq 'End')
	 )
	{
#	  $overlapOK=1;
	  debug(10,"No overlap Actual due to $strandActual $alignmentActual --- ignored...  ")
	}

  #### ckeck now the second contig
  if (
	  ($strand2 eq '+')
	 ){
	if ($start2 < $MAX_CONTIG_END)  {
	  $alignmentNew='Begin';
	  debug(10,"New Overlap + Begin");
	}
	elsif ($end2 > ($lenNew-$MAX_CONTIG_END)) {
	  $alignmentNew='End';
	  debug(10,"New Overlap + End");
	}
	else {
	  debug(10,"Deeep shit New + ");
	  # $overlapOK=0
	}
  }
  else {
	### overlap is -
	if ($end2 < $MAX_CONTIG_END)  {
	  debug(10,"New Overlap -Begin");
	  $alignmentNew='Begin';
	}
	elsif ($start2 > ($lenNew-$MAX_CONTIG_END)) {
	  debug(10,"New Overlap - End")	;

	  $alignmentNew='End';
	}
	else {
	  debug(10,"Deeep shit New - ");
	  #$overlapOK=0
	}
  }
  if (
	  ($strandNew==1 && $alignmentNew eq 'End')
	  or
	  ($strandNew==-1 && $alignmentNew eq 'Start')
	 )
	{
	  $overlapOK=0;
	  debug(10,"No overlap -- New ")
	}


  return $overlapOK;

}
sub makeN{
  my $num    = shift;
  my $str='';

  for (1..$num) {
	$str.="N";

  }

  return "$str\n";

}
sub revcomp{
  my $str = shift;

  my $str2= reverse($str);

	$str2 =~ tr/atgcATGC/tacgTACG/;


  return $str2;
}
=head2

        Usage           :

        Arg [1]         : LSF queue for jobs

        Arg [2]         : file of mapped lanes to be summarised

        Arg [3]         : faidx index file for reference

        Example         :

        Description     : Summarises the coverage of the mapped lanes specified in the lanes.fofn file

        Returntype      : none

        Author          : Jacqueline McQuillan E<lt>jm15@sanger.ac.uk<gt>

=cut

sub get_contig_coverage
{

#code for the function goes here

}

=head1 AUTHOR

Team 133, E<lt>team133g@sanger.ac.uk<gt>

=head1 COPYRIGHT AND LICENSE

Copyright (C) 2009 Genome Research Limited. All Rights Reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

=cut

1;
__END__
