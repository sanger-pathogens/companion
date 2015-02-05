#!/bin/bash

chr=$1;
path=$2

if [ -z "$path" ] ; then

	if [ -f "embl/$chr.embl" ] ; then
		act embl/$chr.embl comp/comp.$chr.blast <(cat Res.$chr.contigs.gff Res.$chr.fna)
	else
		act Reference/$chr comp/comp.$chr.blast <(cat Res.$chr.contigs.gff Res.$chr.fna)
	fi
else

### the "first comp of the top genome must be switch, as subject and query are revcomp"
cat $path/comp/comp.$chr.blast | perl -nle '@ar=split(/\t/);$tmp=$ar[6];$ar[6]=$ar[8];$ar[8]=$tmp;$tmp=$ar[7];$ar[7]=$ar[9];$ar[9]=$tmp;print join("\t",@ar)' > $path/comp/REV.comp.$chr.blast

 if [ -f "embl/$chr.embl" ] ; then
echo "act <(cat $path/Res.$chr.contigs.gff  $path/Res.$chr.fna) $path/comp/REV.comp.$chr.blast embl/$chr.embl comp/comp.$chr.blast <(cat Res.$chr.contigs.gff Res.$chr.fna)"
  act <(cat $path/Res.$chr.contigs.gff  $path/Res.$chr.fna) $path/comp/REV.comp.$chr.blast embl/$chr.embl comp/comp.$chr.blast <(cat Res.$chr.contigs.gff Res.$chr.fna)
 else
echo "act <(cat $path/Res.$chr.contigs.gff  $path/Res.$chr.fna) $path/comp/REV.comp.$chr.blast Reference/$chr comp/comp.$chr.blast <(cat Res.$chr.contigs.gff Res.$chr.fna)"
 act <(cat $path/Res.$chr.contigs.gff  $path/Res.$chr.fna) $path/comp/REV.comp.$chr.blast Reference/$chr comp/comp.$chr.blast <(cat Res.$chr.contigs.gff Res.$chr.fna)
 fi
fi