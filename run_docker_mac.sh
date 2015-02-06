#!/bin/sh
export ROOTDIR=`pwd`
./nextflow run annot.nf -resume -with-docker satta/annot-nf
