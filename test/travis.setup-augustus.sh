#!/bin/bash
set -ev
wget -q http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.0.3.tar.gz
tar -xzvf augustus*.tar.gz
rm -rf augustus*.tar.gz
mv augustus* augustus
rm -rf augustus/docs