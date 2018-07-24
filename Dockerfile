#
#  From this base-image / starting-point
#
FROM debian:testing

#
#  Authorship
#
MAINTAINER ss34@sanger.ac.uk

#
# Pull in packages from testing
#
RUN apt-get update -qq

#
# Install dependencies
# we need blast2 for ABACAS2 as it does not use BLAST+ yet
#
RUN apt-get install build-essential hmmer lua5.1 ncbi-blast+ blast2 snap \
                    unzip mummer infernal exonerate mafft fasttree \
                    circos libsvg-perl libgd-svg-perl python-setuptools \
                    libc6-i386 lib32stdc++6 lib32gcc1 netcat genometools \
                    last-align libboost-iostreams-dev libgslcblas0 libgsl-dev \
                    libcolamd2 liblpsolve55-dev libstdc++6 aragorn tantan \
                    libstorable-perl libbio-perl-perl libsqlite3-dev \
                    --yes
RUN ln -fs /usr/bin/fasttree /usr/bin/FastTree

#
# Install AUGUSTUS
#

# additional dependencies for 3.3.1
RUN apt-get install libboost-graph-dev libsuitesparse-dev zlib1g-dev libbam-dev --yes # libhts-dev

# Augustus 3.3.1 includes bam2hints and bam2wig which have dependencies
# bam2hints dependencies
RUN apt-get install bamtools libbamtools-dev --yes
# bam2wig dependencies
# htslib, bcftools, samtools & tabix must be downloaded & built from github, rather than installed from debian repo
RUN apt-get install autoconf zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev libcrypto++-dev libssl-dev libncurses5-dev --yes 
ADD https://github.com/samtools/htslib/archive/develop.zip /opt/htslib-develop.zip
RUN cd /opt && unzip htslib-develop.zip && cd /opt/htslib-develop && \
    autoheader && autoconf && ./configure && \
    make && make install && \
    rm -rf /opt/htslib-develop*
ADD https://github.com/samtools/bcftools/archive/develop.zip  /opt/bcftools-develop.zip
RUN cd /opt && unzip bcftools-develop.zip && cd /opt/bcftools-develop && \
    autoheader && autoconf && ./configure && \
    make && make install && \
    rm -rf /opt/bcftools-develop*
ADD https://github.com/samtools/samtools/archive/develop.zip  /opt/samtools-develop.zip
RUN cd /opt && unzip samtools-develop.zip && cd /opt/samtools-develop && \
    autoheader && autoconf -Wno-syntax && ./configure && \
    make && make install && \
    rm -rf /opt/samtools-develop*
ADD https://github.com/samtools/tabix/archive/master.zip /opt/tabix-master.zip
RUN cd /opt && unzip tabix-master.zip && cd /opt/tabix-master && \
    make && \
    rm -rf /opt/tabix-master*

# augustus compile looks in non-existence /root/tools for things that just got installed in /usr/local/include
# temporary (yeah, right...) workaround:
RUN ln -s /usr/local/include /root/tools

ADD http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.3.1.tar.gz /opt/augustus-3.3.1.tar.gz
RUN cd /opt && \
    tar -xzvf augustus* && \
    rm -rf *.tar.gz && \
    mv augustus* augustus && \
    rm -rf augustus/docs && \
    cd augustus && \
    make -j2

#
# Install GenomeTools (most recent git master)
#
ADD https://github.com/genometools/genometools/archive/master.zip /opt/genometools-master.zip
RUN cd /opt && \
    unzip genometools-master.zip && \
    cd /opt/genometools-master && \
    make -j3 cairo=no curses=no prefix=/usr && \
    make -j3 cairo=no curses=no prefix=/usr install && \
    cd gtpython && \
    python setup.py install && \
    cd / && \
    rm -rf /opt/genometools-master*

#
# Install and configure OrthoMCL
#
ADD http://www.orthomcl.org/common/downloads/software/unsupported/v1.4/ORTHOMCL_V1.4_mcl-02-063.tar /opt/omcl.tar
RUN cd /opt && \
    tar -xvf omcl.tar && \
    tar -xzvf mcl-02-063.tar.gz && \
    rm -f omcl.tar mcl-02-063.tar.gz && \
    cd /opt/mcl-* && \
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
# Install Gblocks
#
ADD http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux64_0.91b.tar.Z /opt/gblocks64.tar.Z
RUN cd /opt && \
    tar -xzvf gblocks64.tar.Z && \
    rm -rf gblocks64.tar.Z && \
    cp Gblocks_0.91b/Gblocks /usr/bin/Gblocks && \
    chmod 755 /usr/bin/Gblocks

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

ENV AUGUSTUS_CONFIG_PATH /opt/augustus/config
ENV RATT_HOME /opt/RATT
ENV GT_RETAINIDS yes
ENV PERL5LIB /opt/ORTHOMCLV1.4/:/opt/RATT/:/opt/ABACAS2/:$PERL5LIB
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/augustus/bin:/opt/augustus/scripts:/opt/ORTHOMCLV1.4:/opt/RATT:/opt/ABACAS2:$PATH
