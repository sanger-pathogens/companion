#
#  From this base-image / starting-point
#
FROM debian:wheezy

#
#  Authorship
#
MAINTAINER ss34@sanger.ac.uk

#
# Pull in packages from testing
#
RUN echo "deb http://http.debian.net/debian testing main" > /etc/apt/sources.list

#
# Update apt
#
RUN apt-get update -q -q

#
# Install dependencies from Debian
#
RUN apt-get install build-essential hmmer lua5.1 blast2 snap \
                    liblua5.1-0 libcairo2 zlib1g libbz2-1.0 libexpat1 libpth20 \
                    libncurses5 libsqlite3-0 libpango-1.0-0 \
                    libpangocairo-1.0-0 libtre5 python-ctypes liblua5.1-0-dev \
                    lua-md5-dev lua-filesystem-dev lua-lpeg-dev libcairo2-dev \
                    zlib1g-dev libbz2-dev libexpat1-dev libncurses5-dev \
                    libsqlite3-dev libbam-dev libpango1.0-dev libtre-dev \
                    python unzip libstorable-perl cpanminus mummer \
                    infernal exonerate \
                    --yes --force-yes

#
# Install AUGUSTUS (binaries)
#
ADD http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.0.3.tar.gz /opt/augustus-3.0.3.tar.gz
RUN cd /opt && \
    tar -xzvf augustus* && \
    rm -rf *.tar.gz && \
    mv augustus* augustus

#
# Install GenomeTools (most recent git master)
#
ADD https://github.com/genometools/genometools/archive/master.zip /opt/genometools-master.zip
RUN cd /opt && \
    unzip genometools-master.zip && \
    cd /opt/genometools-master && \
    make -j3 && \
    make -j3 install && \
    cd / && \
    rm -rf /opt/genometools-master*

#
# Install ARAGORN
# (this is done from source instead of Debian b/c Debian-built version hangs)
#
ADD http://mbio-serv2.mbioekol.lu.se/ARAGORN/Downloads/aragorn1.2.36.tgz /opt/aragorn.tgz
RUN cd /opt && \
    tar -xvf aragorn.tgz && \
    mv aragorn1* aragorn && \
    cd aragorn && \
    gcc -O2 -o aragorn aragorn*.c

#
# Add Perl deps (needed for OrthoMCL)
#
RUN cpanm --force Carp Bio::SearchIO List::Util Getopt::Long && \
    rm -rf /root/.cpanm/work/

#
# Install and configure OrthoMCL
#
ADD http://www.orthomcl.org/common/downloads/software/unsupported/v1.4/ORTHOMCL_V1.4_mcl-02-063.tar /opt/omcl.tar
RUN cd /opt && \
    tar -xvf omcl.tar && \
    tar -xzvf mcl-02-063.tar.gz && \
    rm -f omcl.tar mcl-02-063.tar.gz
RUN cd /opt/mcl-* && \
    ./configure && \
    make -j3 && \
    make install && \
    cd / && \
    rm -rf /opt/mcl*
RUN sed -i 's/our .PATH_TO_ORTHOMCL.*=.*/our $PATH_TO_ORTHOMCL = ".\/";/' /opt/ORTHOMCLV1.4/orthomcl_module.pm && \
    sed -i 's/our .BLASTALL.*=.*/our $BLASTALL = "\/usr\/bin\/blastall";/' /opt/ORTHOMCLV1.4/orthomcl_module.pm && \
    sed -i 's/our .FORMATDB.*=.*/our $FORMATDB = "\/usr\/bin\/formatdb";/' /opt/ORTHOMCLV1.4/orthomcl_module.pm && \
    sed -i 's/our .MCL.*=.*/our $MCL = "\/usr\/local\/bin\/mcl";/' /opt/ORTHOMCLV1.4/orthomcl_module.pm

#
# get GO OBO file
#
ADD http://geneontology.org/ontology/go.obo /opt/go.obo

#
# get Pfam pHMMs
#
RUN mkdir -p /opt/pfam
ADD http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz /opt/pfam/Pfam-A.hmm.gz
RUN cd /opt/pfam && \
    gunzip Pfam-A.hmm.gz && \
    hmmpress Pfam-A.hmm && \
    rm -f Pfam-A.hmm

#
# copy data dir
#
RUN mkdir -p /opt/data
ADD ./data /opt/data

#
# install RATT (keep up to date from build directory)
#
ADD ./RATT /opt/RATT

#
# install ABACAS (keep up to date from build directory)
#
ADD ./ABACAS2 /opt/ABACAS2

#
# clean up dev stuff (not strictly necessary)
#
RUN apt-get purge build-essential --yes --force-yes
RUN apt-get remove liblua5.1-0-dev lua-md5-dev lua-filesystem-dev lua-lpeg-dev \
                   libcairo2-dev zlib1g-dev libbz2-dev libexpat1-dev \
                   libncurses5-dev libsqlite3-dev libbam-dev libpango1.0-dev \
                   libtre-dev --yes --force-yes

ENV AUGUSTUS_CONFIG_PATH /opt/augustus/config
ENV RATT_HOME /opt/RATT
ENV GT_RETAINIDS yes
ENV PERL5LIB /opt/ORTHOMCLV1.4/:/opt/RATT/:/opt/ABACAS2/:$PERL5LIB
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/augustus/bin:/opt/augustus/scripts:/opt/ORTHOMCLV1.4:/opt/RATT:/opt/ABACAS2:/opt/aragorn:$PATH
