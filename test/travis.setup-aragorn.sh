#!/bin/bash
set -ev
wget http://mbio-serv2.mbioekol.lu.se/ARAGORN/Downloads/aragorn1.2.36.tgz
tar -xvf aragorn*.tgz
rm aragorn*.tgz
mv aragorn* aragorn
cd aragorn
gcc -O2 -o aragorn aragorn*.c