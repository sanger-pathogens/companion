#!/bin/bash
set -ev
wget -q http://mbio-serv2.mbioekol.lu.se/ARAGORN/Downloads/aragorn1.2.36.tgz
tar -xvf aragorn*.tgz
rm aragorn*.tgz
mv aragorn* aragorn
cd aragorn
gcc -Os -o aragorn aragorn*.c
