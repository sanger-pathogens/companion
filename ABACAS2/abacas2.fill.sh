#!/bin/bash

reference=$1
contig=$2
ABA_MIN_LENGTH=$3
ABA_MIN_IDENTITY=$4
doBlast=$5
MemBig=$6

if [ -z "$contig" ] ; then
	echo "

                 *** Abacas II. ***       For any distrubation with this program, please don't blame Sammy!

usage:
~tdo/Bin/abacas2.sh <reference> <Contig to order> optinal: <min aligment length> <Identity cutoff> <doblast: 0/1> <Bsub RAM>

reference:           Fasta (or multi-fasta) against which the contigs should be orders
contig:              Contigs or query that should be ordered
Min aligment length: Threshold for the length, when an alignment lenght is significant. (default 200)
Identity cutoff:     Threshold for identity to place contigs. (default 95)
Do Blast:            Does a blast for the act. (default 0)
Bsub RAM:            KB used for the abacas2 run in the second stage. (default: 6000) Attention, if run of farm2, more than 12000 shouldn't be used, as it should go to hugemem.

Further parameters:
ABA_CHECK_OVERLAP=0; export ABA_CHECK_OVERLAP # this will not try to overlap contigs, quite quicker 10teims
ABA_splitContigs=1; export ABA_splitContigs # this parameter will split contigs. This is good to split the orign, and to find rearrangement. A split contigs has the suffix _i (i the part)
ABA_WORD_SIZE  # sets the word size. This is critical for speed issues in nucmer. default is 20
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
   MemBig=6000
fi

if [ "$ABA_MIN_IDENTITY" -gt "99" ] ; then
	echo "Your identity might be too high $ABA_MIN_IDENTITY > 99 "
	exit ;
fi
tmp=$$
sed 's/|/_/g' $reference > Ref.$tmp
reference=Ref.$tmp
ln -s $contig Contigs.$tmp
contig=Contigs.$tmp


ABA_WORD_SIZE="50 --maxmatch "; export ABA_WORD_SIZE

export ABA_MIN_LENGTH ABA_MIN_IDENTITY contig reference

tmp=$$

echo "~tdo/Bin/abacasII/run.AbacasII.comparison.sh $reference $contig; perl ~tdo/Bin/abacas2.getFillSeq.pl $contig.coords $reference >> $contig"
bsub -J"stage0.$tmp" -R "select[type==X86_64 && mem > 8000] rusage[mem=8000]" -M8000000 -o abacas2_0.o -e abacas2_0.e "~tdo/Bin/abacasII/run.AbacasII.comparison.sh $reference $contig; perl ~tdo/Bin/abacas2.getFillSeq.pl $contig.coords $reference >> $contig"


bsub -J"stage1.$tmp"  -w"stage0.$tmp" -R "select[type==X86_64 && mem > 8000] rusage[mem=8000]" -M8000000 -o abacas2.o -e abacas2.e "~tdo/Bin/abacasII/run.AbacasII.comparison.sh $reference $contig"


#do ordering
for x in `grep '>' $reference | perl -nle '/>(\S+)/;print $1' ` ;
  do
#  echo " bsub -w\"stage1.$tmp\" -J\"stage2.$tmp\" -R \"select[type==X86_64 && mem > $MemBig] rusage[mem=$MemBig]\" -M\"$MemBig\"000 -o output2.$x.o -e output2.$x.e \"~tdo/Bin/goTilingGraph.pl $x.coords  $contig Res\";"
  bsub -w"stage1.$tmp" -J"stage2.$tmp" -R "select[type==X86_64 && mem > $MemBig] rusage[mem=$MemBig]" -M"$MemBig"000 -o output2.$x.o -e output2.$x.e "~tdo/Bin/goTilingGraph.pl $x.coords  $contig Res";
done


#do bin
bsub -w "ended(stage2.$tmp)" -J"stage3.$tmp" -R "select[type==X86_64 && mem > 6000] rusage[mem=6000]" -M6000000 -o output2.$x.o -e output2.$x.e "~tdo/Bin/abacas2.bin.sh $contig Res.abacasBin.fna && grep -v '>'  Res.abacasBin.fna | awk 'BEGIN {print \">Bin.union\"} {print}' > Res.abacasBin.oneSeq.fna_ && cat Res*fna > tmp && /software/pathogen/projects/protocols/bin/fasta2singleLine.pl tmp tmp2  && /nfs/users/nfs_j/jit/repository/pathogen/user/jit/converter/fasta2multipleLine_v2.pl tmp2 Genome.abacas.fasta 80  &> /dev/null && mv  Res.abacasBin.fna_ Res.abacasBin.fna && rm tmp tmp2";


if [ "$doBlast" == "1" ] ; then
	bsub -w "ended(stage2.$tmp)" -J"stage3.$tmp.blast"	~tdo/Bin/abacas2.doblast.sh $reference Res
fi

