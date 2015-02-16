#!/bin/bash
set -ev
wget http://mbio-serv2.mbioekol.lu.se/ARAGORN/Downloads/aragorn1.2.36.tgz
tar -xvf aragorn1.2.36.tgz
mv aragorn1* aragorn
cd aragorn
gcc -O2 -o aragorn aragorn*.c