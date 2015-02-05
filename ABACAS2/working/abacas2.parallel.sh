#!/bin/bash

reference=$1
contig=$2
ABA_MIN_LENGTH=$3
ABA_MIN_IDENTITY=$4
doBlast=$5
MemBig=$6

PERL5LIB=$PERL5LIB:/nfs/pathogen003/tdo/Tools/ABACAS2/working
export PERL5LIB

if [ ! -f "$contig" ] ; then
	echo "

                 *** Abacas II. ***       For any distrubation with this program, please don't blame Sammy!

usage:
abacas2.sh <reference> <Contig to order> optinal: <min aligment length> <Identity cutoff> <doblast: 0/1> <Bsub RAM>

reference:           Fasta (or multi-fasta) against which the contigs should be orders
contig:              Contigs or query that should be ordered
Min aligment length: Threshold for the length, when an alignment lenght is significant. (default 200)
Identity cutoff:     Threshold for identity to place contigs. (default 95)
Do Blast:            Does a blast for the act. (default 1)
Bsub RAM:            KB used for the abacas2 run in the second stage. (default: 6000) Attention, if run of farm2, more than 12000 shouldn't be used, as it should go to hugemem.

Further parameters:
ABA_CHECK_OVERLAP=1; export ABA_CHECK_OVERLAP # this will try to overlap contigs, quite slow 10 times
ABA_splitContigs=1; export ABA_splitContigs # this parameter will split contigs. This is good to split the orign, and to find rearrangement. A split contigs has the suffix _i (i the part)
ABA_WORD_SIZE  # sets the word size. This is critical for speed issues in nucmer. default is 20

Advanced (for pipeline usage):
ABACAS_WAIT_START    Will give the first job of abacas a name
ABACAS_LAST_START    Will give the last job a name.
"

exit;
fi

if [ -z "$ABA_MIN_LENGTH" ] ; then
	ABA_MIN_LENGTH=200;
fi
if [ -z "$ABA_MIN_IDENTITY" ] ; then
	ABA_MIN_IDENTITY=95;
fi
if [ -z "$MemBig" ]; then
   MemBig=4500
fi
if [ -z "$ABA_CHECK_OVERLAP" ] ; then
        ABA_CHECK_OVERLAP=0;
        export ABA_CHECK_OVERLAP
fi
if [ -z "$doBlast" ]; then
   doBlast=1
fi

### for pipelining
tmp=$$
if [ ! -z $ABACAS_WAIT_START ] ; then
   ABACAS_WAIT_START=" -w $ABACAS_WAIT_START "
fi
if [ -z $ABACAS_LAST_START ] ; then
   ABACAS_LAST_START="stage3.$tmp"
fi


if [ "$ABA_MIN_IDENTITY" -gt "99" ] ; then
	echo "Your identity might be too high $ABA_MIN_IDENTITY > 99 "
	exit ;
fi

sed 's/|/_/g' $reference > Ref.$tmp
reference=Ref.$tmp
ln -s $contig Contigs.$tmp
contig=Contigs.$tmp

export ABA_MIN_LENGTH ABA_MIN_IDENTITY contig reference


tmp=$$
bsub $ABACAS_WAIT_START -J"stage1.$tmp" \
  -R "select[type==X86_64 && mem > 4000] rusage[mem=4000]" -M4000 \
  -o abacas2.o -e abacas2.e \
  "abacas2.runComparison.sh $reference $contig"


#do ordering, construct job dependencies
i=0
BSUB_W=""
for x in `grep '>' $reference | perl -nle '/>(\S+)/;print $1' `;
do
  ((i++))
  bsub -w"stage1.$tmp" -J"stage2.$tmp.$i" \
    -R "select[type==X86_64 && mem > $MemBig] rusage[mem=$MemBig]" \
    -M"$MemBig" -o output2.$x.o -e output2.$x.e \
    "abacas2.doTilingGraph.pl $x.coords  $contig Res";
  BSUB_W="${BSUB_W}ended(stage2.${tmp}.${i}) && "
done
# cut off last '&&'
BSUB_W=`echo "$BSUB_W" | sed 's/&& $//'`
echo "bsub: $BSUB_W"

if [ "$doBlast" == "1" ] ; then
  bsub -w "$BSUB_W" -J"stage4.$tmp.blast" -o blast.started.o \
    -e blast.started.e  -R "select[type==X86_64 && mem > 500] rusage[mem=500]" \
    -M500  "abacas2.doblast.sh $reference Res"
fi

#do bin
bsub -w "$BSUB_W" -J"$ABACAS_LAST_START" \
  -R "select[type==X86_64 && mem > 6000] rusage[mem=6000]" -M6000 \
  -o output2.bin.o -e output2.bin.e \
  "~tdo/Bin/abacas2.bin.sh $contig Res.abacasBin.fna && grep -v '>'  Res.abacasBin.fna | awk 'BEGIN {print \">Bin.union\"} {print}' > Res.abacasBin.oneSeq.fna_ && cat Res*fna > Genome.abacas.fasta && bam.correctLineLength.sh Genome.abacas.fasta ";


