#!/bin/bash
set -ev
wget https://github.com/genometools/genometools/archive/master.zip
unzip master.zip
cd genometools-master
make -j3 cairo=no curses=no
make -j3 cairo=no curses=no install
cd ..
rm -rf /opt/genometools-master*