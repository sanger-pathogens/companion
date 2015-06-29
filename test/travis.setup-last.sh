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

rm -f tantan-13.zip
wget -q http://cbrc3.cbrc.jp/~martin/tantan/tantan-13.zip
unzip tantan-13.zip
cd tantan-13
make
make install
cd ..
rm -rf tantan-13
rm -f tantan-13.zip