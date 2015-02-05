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
RUN apt-get install build-essential aragorn hmmer lua5.1 ncbi-blast+ snap \
                    liblua5.1-0 libcairo2 zlib1g libbz2-1.0 libexpat1 libpth20 \
                    libncurses5 libsqlite3-0 libpango-1.0-0 \
                    libpangocairo-1.0-0 libtre5 python-ctypes liblua5.1-0-dev \
                    lua-md5-dev lua-filesystem-dev lua-lpeg-dev libcairo2-dev \
                    zlib1g-dev libbz2-dev libexpat1-dev libncurses5-dev \
                    libsqlite3-dev libbam-dev libpango1.0-dev libtre-dev \
                    python unzip bioperl libstorable-perl mummer \
                    --yes --force-yes

#
# Install AUGUSTUS
#
ADD http://augustus.gobics.de/binaries/augustus.2.5.5.tar.gz /opt/augustus.2.5.5.tar.gz
RUN cd /opt && \
    tar -xzvf augustus* && \
    rm -rf *.tar.gz && \
    mv augustus* augustus && \
    cd augustus && \
    make

#
# Install GenomeTools (most recent git master)
#
ADD https://github.com/genometools/genometools/archive/master.zip /opt/genometools-master.zip
RUN cd /opt && \
    unzip genometools-master.zip && \
    cd /opt/genometools-master && \
    make -j3 && \
    make -j3 install && \
    make cleanup && \
    cd / && \
    rm -rf /opt/genometools-master*

#
# Install OrthoMCL
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

#
# get GO OBO file
#
ADD http://geneontology.org/ontology/go.obo /opt/go.obo

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

# TODO: add ref species models etc.
# ADD fgram_base /opt/augustus/config/species/

ENV AUGUSTUS_CONFIG_PATH /opt/augustus/config
ENV RATT_HOME /opt/RATT
ENV GT_RETAINIDS yes
ENV PERL5LIB /opt/ORTHOMCLV1.4/:$PERL5LIB
ENV PATH_TO_ORTHOMCL /opt/ORTHOMCLV1.4/
ENV PATH /usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/augustus/bin:/opt/augustus/scripts:/opt/ORTHOMCLV1.4:/opt/RATT:/opt/ABACAS2
