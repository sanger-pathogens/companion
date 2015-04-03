#!/bin/bash
set -ev
rm -f master.zip
wget -q https://github.com/sanger-pathogens/gff3toembl/archive/master.zip
unzip master.zip
cd gff3toembl-master
python setup.py install
cd ..
rm -rf /opt/gff3toembl-master*
rm -f master.zip