#!/bin/bash
set -ev
wget https://github.com/sanger-pathogens/gff3toembl/archive/master.zip
unzip master.zip
cd gff3toembl-master
python setup.py install
cd ..
rm -rf /opt/gff3toembl-master*
rm master.zip