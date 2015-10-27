#!/bin/bash
set -ev
wget -q http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.2.tar.gz
tar -xzvf augustus*.tar.gz
rm -rf augustus*.tar.gz
mv augustus* augustus
rm -rf augustus/docs
