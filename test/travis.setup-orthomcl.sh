#!/bin/bash
set -ev
wget -q http://www.orthomcl.org/common/downloads/software/unsupported/v1.4/ORTHOMCL_V1.4_mcl-02-063.tar
tar -xvf ORTHOMCL_V1.4_mcl-02-063.tar
tar -xzvf mcl-02-063.tar.gz
rm -f ORTHOMCL_V1.4_mcl-02-063.tar mcl-02-063.tar.gz
cd mcl-*
./configure
make -j3
sudo make install
cd ..
sed -i 's/our .PATH_TO_ORTHOMCL.*=.*/our $PATH_TO_ORTHOMCL = ".\/";/' ./ORTHOMCLV1.4/orthomcl_module.pm
sed -i 's/our .BLASTALL.*=.*/our $BLASTALL = "\/usr\/bin\/blastall";/' ./ORTHOMCLV1.4/orthomcl_module.pm
sed -i 's/our .FORMATDB.*=.*/our $FORMATDB = "\/usr\/bin\/formatdb";/' ./ORTHOMCLV1.4/orthomcl_module.pm
sed -i 's/our .MCL.*=.*/our $MCL = "\/usr\/local\/bin\/mcl";/' ./ORTHOMCLV1.4/orthomcl_module.pm