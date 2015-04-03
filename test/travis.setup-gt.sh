#!/bin/bash
set -ev
rm -f master.zip
wget -q https://github.com/genometools/genometools/archive/master.zip
unzip master.zip
cd genometools-master
make -j3 cairo=no curses=no prefix=/usr
make -j3 cairo=no curses=no prefix=/usr install
cd gtpython
python setup.py install
cd ..
rm -rf /opt/genometools-master*
rm -f master.zip