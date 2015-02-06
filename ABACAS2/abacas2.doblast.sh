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
SeperateSequences.pl  ../REF.$tmp
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
echo  " blastall -p blastn -F F -W 15 -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast"
   blastall -p blastn -W 25 -F T -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast
 else

echo  " blastall -p blastn  $Other -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast"
   blastall -p blastn  $Other -m 8 -e 1e-20 -d Reference/$nameRes -i $pre.$nameRes.fna -o comp/comp.$nameRes.blast
fi
done
