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
ABA_CHECK_OVERLAP=1; export ABA_CHECK_OVERLAP # this will not try to overlap contigs, quite quicker 10teims
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
if [ -z "$ABA_CHECK_OVERLAP" ] ; then
        ABA_CHECK_OVERLAP=0;
	export ABA_CHECK_OVERLAP
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

export ABA_MIN_LENGTH ABA_MIN_IDENTITY contig reference

tmp=$$


abacas2.runComparison.sh $reference $contig

for x in `grep '>' $reference | perl -nle '/>(\S+)/;print $1' ` ; do
	abacas2.doTilingGraph.pl $x.coords  $contig Res
done

abacas2.bin.sh $contig Res.abacasBin.fna && grep -v '>'  Res.abacasBin.fna | awk 'BEGIN {print ">Bin.union"} {print}' > Res.abacasBin.oneSeq.fna_ && cat Res*fna > Genome.abacas.fasta && bam.correctLineLength.sh Genome.abacas.fasta  &> /dev/null
if [ -f "Res.abacasBin.fna_" ] ; then
    mv Res.abacasBin.fna_ Res.abacasBin.fna
fi

#~tdo/Bin/abacas2.doblast.sh $reference Res
#ef=$reference
#pre=Res

#if [ -z "$pre" ] ; then
#   echo "usage <reference> Res"

#fi

#tmp=$$

#ln -s $ref REF.$tmp

#mkdir Reference
#cd Reference
#SeperateSequences.pl  ../REF.$tmp
#cd ..
#mkdir comp

#count=0
#for nameRes in `grep '>' $ref | perl -nle 's/\|/_/g;/>(\S+)/; print $1'` ; do
#	let count++;
#	if [ $count -gt 200 ] ; then
#		echo "too many contigs to bsub!\n";
#		exit
#	fi
#	formatdb -p F -i Reference/$nameRes;

#	megablast -F T -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast

#if [ -z "$Other" ] ; then
#echo  " bsub -o blast.$nameRes.o -e blast.$nameRes.e -R \"select[type==X86_64 && mem > 6000] rusage[mem=6000]\" -M6000000  blastall -p blastn -F F -W 15 -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast"
#megablast -F T -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast
 #else

 #blastall -p blastn  $Other -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast

#done



