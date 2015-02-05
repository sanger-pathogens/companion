#! /usr/bin/perl -w
#
# File: goTilingGraph.pl
# Time-stamp: <22-Dec-2014 09:18:32 tdo>
# $Id: $
#
# Copyright (C) 2009 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
#

use strict;

use lib '/nfs/pathogen003/tdo/Tools/ABACAS2/working';
use TilingGraph;

if (scalar(@ARGV) < 3) {
  die "Usage: coords sequence Resultname\n";
  
}
TilingGraph::startTiling(shift,shift,shift);




