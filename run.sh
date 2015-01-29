#!/bin/sh
NXF_WORK=/lustre/scratch109/sanger/ss34/nxf-tmp NXF_TEMP=/lustre/scratch109/sanger/ss34/nxf-tmp ../nextflow  -C annot.conf run annot.nf -resume -with-trace
