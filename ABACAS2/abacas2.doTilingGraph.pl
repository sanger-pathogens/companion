#! /usr/bin/perl -w
#
# File: goTilingGraph.pl
# Time-stamp: <01-Oct-2009 11:41:27 tdo>
# $Id: $
#
# Copyright (C) 2009 by Pathogene Group, Sanger Center
#
# Author: Thomas Dan Otto
#
# Description:
#

use strict;

use TilingGraph;

if (scalar(@ARGV) < 3) {
  die "Usage: coords sequence Resultname\n";

}
TilingGraph::startTiling(shift,shift,shift);




