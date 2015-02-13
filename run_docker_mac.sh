#!/bin/sh
export ROOTDIR=`pwd`
nextflow run -c loc_docker.config -c params_default.config annot.nf -resume -with-docker satta/annot-nf
