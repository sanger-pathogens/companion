#!/bin/bash
set -ev
rm -f last-581.zip
wget -q http://last.cbrc.jp/last-581.zip
unzip last-581.zip
cd last-581
make
make install
cd ..
rm -rf last-581
rm -f last-581.zip
