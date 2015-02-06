#!/bin/bash

contig=$1
resultName=$2
pre=$3

if [ -z "$pre" ] ; then
	pre=Res
fi

#echo "grep Contig $pre.*.contigs.gff | perl -nle '/label=\"\S*contig=(\S+)\";note/;print $1' > List.mappingcontig.fofn"

grep Contig $pre.*.contigs.gff | perl -nle '/label=\"\S*contig=(\S+)\";note/;print $1' > List.mappingcontig.fofn

Assembly.deleteContigs.pl List.mappingcontig.fofn $contig $resultName
#rm  List.mappingcontig.fofn
