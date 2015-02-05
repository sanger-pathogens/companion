#! /bin/bash
#
# File: bam.correctLineLength.pl
# Time-stamp: <10-Aug-2010 14:34:53 tdo>
# $Id: $
#
# Copyright (C) 2010 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
# 



name=$1

if [ -z "$name" ] ; then
	echo "Please give the (multi-) fasta to be corrected...";
exit;
fi

tmp=$$

PERL5LIB=$PERL5LIB:/software/pathogen/projects/protocols/lib/perl5
export PERL5LIB

cat $name | perl ~tdo/Bin/fasta2singleLine.tdo.pl > $name.$tmp
#/nfs/users/nfs_j/jit/repository/pathogen/user/jit/converter//
fasta2multipleLine_v2.pl $name.$tmp $name 80
rm $name.$tmp
