#!/bin/sh
export ROOTDIR=`pwd`
NXF_WORK=/lustre/scratch109/sanger/ss34/nxf-tmp NXF_TEMP=/lustre/scratch109/sanger/ss34/nxf-tmp ../nextflow run annot.nf
