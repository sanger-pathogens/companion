#!/bin/bash

ref=$1
pre=$2
Other=$3

if [ -z "$pre" ] ; then
   echo "usage <reference> Res"

fi

tmp=$$

ln -s $ref REF.$tmp

mkdir Reference
cd Reference
/nfs/users/nfs_t/tdo/Bin/SeperateSequences.pl  ../REF.$tmp
cd ..
mkdir comp
 
count=0
for nameRes in `grep '>' $ref | perl -nle 's/\|/_/g;/>(\S+)/; print $1'` ; do
   let count++;
   if [ $count -gt 200 ] ; then
	  echo "too many contigs to bsub!\n";
	  exit
   fi
   formatdb -p F -i Reference/$nameRes;

if [ -z "$Other" ] ; then
echo  " bsub -o blast.$nameRes.o -e blast.$nameRes.e -R \"select[type==X86_64 && mem > 6000] rusage[mem=6000]\" -M6000000  blastall -p blastn -F F -W 15 -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast"
   bsub -o blast.$nameRes.o -e blast.$nameRes.e -R "select[type==X86_64 && mem > 6000] rusage[mem=6000]" -M6000  blastall -p blastn -W 25 -F T -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast
 else

echo  " bsub -o blast.$nameRes.o -e blast.$nameRes.e -R \"select[type==X86_64 && mem > 6000] rusage[mem=6000]\" -M6000000  blastall -p blastn  $Other -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast"
   bsub -o blast.$nameRes.o -e blast.$nameRes.e -R "select[type==X86_64 && mem > 6000] rusage[mem=6000]" -M6000  blastall -p blastn  $Other -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast
fi
done
