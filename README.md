# Companion
A portable, scalable eukaryotic genome annotation pipeline implemented in Nextflow.

[![Build Status](https://travis-ci.org/sanger-pathogens/companion.svg?branch=master)](https://travis-ci.org/sanger-pathogens/companion)  
[![License: ISC](https://img.shields.io/badge/License-ISC-brightgreen.svg)](https://github.com/sanger-pathogens/companion/blob/master/LICENSE)  
[![status](https://img.shields.io/badge/NAR-10.1093%2Fnar.gkw292-brightgreen.svg)](https://doi.org/10.1093/nar/gkw292)


## Content

## Introduction
This software is a comprehensive computational pipeline for the annotation of eukaryotic genomes (like protozoan parasites). It performs the following tasks:

  - Fast generation of pseudomolecules from scaffolds by ordering and orientating against a reference
  - Accurate transfer of highly conserved gene models from the reference
  - _De novo_ gene finding as a complement to the gene transfer
  - Non-coding RNA detection (tRNA, rRNA, sn(o)RNA, ...)
  - Pseudogene detection
  - Functional annotation (GO, products, ...)
    - ...by transferring reference annotations to the target genome
    - ...by inferring GO terms and products from Pfam pHMM matches
  - Consistent gene ID assignment
  - Preparation of validated GFF3, GAF and EMBL output files for jump-starting manual curation and quick turnaround time to submission

It supports parallelized execution on a single machine as well as on large cluster platforms (LSF, SGE, ...).

## Quick start

This should get you up & running on an Ubuntu system, but please read the full documentation before before doing any work "for real".

### 1. Install dependencies

Execute these commands as root, e.g. using `sudo`
```
apt-get install default-jre
curl -fsSL get.nextflow.io | bash
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add - && \
   add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable" && \
   apt-get update && \
   apt-cache policy docker-ce && \
   apt-get install --yes docker-ce && \
   systemctl enable docker
```

To enable you to use docker with your normal user account (i.e. without being root or needing to sue sudo),
run the following command, with your username in place of `<username>`.  If more than one user will be using _Companion_,
repeat this command with each of their usernames.
```
usermod -aG docker <username>
```

#### Check

  - `java -version` should say you have Java 1.8 or greater
  - `nextflow info` will print system information if nextflow has been installed successfully
  - `systemctl status docker` will tell you if docker is active (running)
  - `docker info` will print information if docker has been instaleld successfully, and you have permission to use it

### 2. Install _Companion_

Execute these command in the directory you want to keep your _Companion_ work in.  Do this as a normal user, i.e. not as root or using sudo.
Use a name that is meaningful to you in place of `<my-companion-project>`
```
curl -L -o companion-master.zip https://github.com/sanger-pathogens/companion/archive/master.zip && \
   unzip companion-master.zip && \
   mv companion-master <my-companion-project>
docker pull sangerpathogens/companion
```

### 3. Run _Companion_ test job

_Companion_ is distributed with configuration and data (including a few pregenerated reference annotations) for a small test run.  Run
the following command (using the name you chose for your project directory in place of `my-companion-project`).
```
nextflow run my-companion-project -profile docker
```
This will create a directory `my-companion-project/example-output` with the results of the run.

### 4. Configure _Companion_ for your annotation run

The file `params_default.config` configures the pipeline, and will need to be edited for your annotation run.  You will probably need to change at least the following parameters:

*inseq*  Your input FASTA file  (`${baseDir}/example-data/L_donovani.1.fasta` in the example parameter file included wirth the distribution)

*ref_dir* The directory containing your reference genomes (`${baseDir}/example-data/references` in the example file)

*ref_species* The "short name" for your reference species (`LmjF.1` in the example file)

*dist_dir* The directory that will contain the newly created output files (`${baseDir}/example-data-output` in the example file)

*GENOME_PREFIX* Text pattern matching your genome prefix (`LDON` in the example file)

*CHR_PATTERN* Pattern matching your chromosome names (`LDON_(%w+)` in the example file, where `%w+` matches one or more letters or numbers)

*ABACAS_BIN_CHR* Abacas bin chromosome <whut?> (`LDON_0` in the example file)

*EMBL_AUTHORS* etc.; please provide suitable EMBL metadata (dummy values in the example file)

*TAXON_ID* Please provide suitable value for the GAF output (`4711` in the example file)

### 5. Prepare reference annotations

The reference annotations used in the pipeline need to be pre-processed before they can be used.  To add a reference organism, you will need:

- a descriptive name of the organism
- a short abbreviation for the organism
- the genome sequence in a single FASTA file
- a structural gene annotation in GFF3 format (see below for details)
- functional GO annotation in GAF 1.0 format, on the gene level
- a pattern matching chromosome headers, describing how to extract chromosome numbers from them
- an [AUGUSTUS](http://bioinf.uni-greifswald.de/augustus/) model, trained on reference genes

Insert these file names, etc., where `<placeholders>` appear in the steps below:

1. Create a new data directory (i.e. the equivalent of the `example-data` directory included in the distribution)
1. Edit `nextflow.config` (and any config files that are referenced) and change parameters such as
`inseq` and `ref_dir` to your new data directory.
1. Copy the new reference genome (FASTA) into `<new_data_dir>/genomes`
1. Copy GFF3 and GAF files into `<new_data_dir>/genomes`
1. Copy Augustus model files into `data/augustus/species/<species_name>/`
1. Create new directory `<new_data_dir>/references/<short_name>/`
1. Add new section to `<new_data_dir>/references/references-in.json`, using the
short name (same as the directory name in the previous step); in this section add
the names/paths of the files copied (above), a descriptive name, and
a pattern for matching chromosomes in the FASTA files (in this example, <short_name>_<n>, where _n_ in any integer).
```
"<short_name>" : {   "gff"                : "../genomes/<gff3_filename>.gff3",
                     "genome"             : "../genomes/<ref_genome_name>.fasta",
                     "gaf"                : "../genomes/<ref_annot_filename>.gaf",
                     "name"               : "<Descriptive Name of Reference Genome>",
                     "augustus_model"     : "../../data/augustus/species/<species_name>/",
                     "chromosome_pattern" : "<short_name>_(%d+)"
                  }
```
8. Finally, change directory to `<new_data_dir>/references` (you _must_ execute the following command in this directory)
and run `../../bin/update_references.lua`.  This writes the file `<new_data_dir>/references/references.json`.

### 6. Run it!

The following command (using the name you chose for your project directory in place of `my-companion-project`) will
start your annotation run:
```
nextflow run my-companion-project -profile docker
```


## Further technical information

### Dependencies

Companion has the following dependencies:

  - [Java](https://openjdk.java.net/) 8 or later
  - [Nextflow](http://nextflow.io)
  - [Docker](https://www.docker.com/) (if using the Docker image to satisfy dependencies)

#### Java

To check if you have Java installed, and the version, use the command `java -version`.  Note that this will give you a version number
of 1.8 for Java 8, 1.9 for Java 9, etc.

To install Java 8 on an Ubuntu or Debian system, run:
```
apt-get install openjdk-8-jre
```
On Fedora, Centos or Red Hat (etc.) systems:
```
yum install java-1.8.0-openjdk
```

#### Nextflow

To install Nextflow, run:
```
curl -fsSL get.nextflow.io | bash
```
This will create an executable called 'nextflow', which should be moved to a suitable directory, for example:
```
mv nextflow /usr/local/bin/
```
Use the command `which nextflow` to check that it is found in your path.

#### Docker

Docker is required if you intended to use the Docker image, as recommended below, to satisfy the dependencies.

To install Docker, see the installation guide for
[Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/),
[Centos](https://docs.docker.com/install/linux/docker-ce/centos/),
[Debian](https://docs.docker.com/install/linux/docker-ce/debian/) or
[Fedora](https://docs.docker.com/install/linux/docker-ce/fedora/).

Users running Companion with Docker will need to be added to the `docker` group (unix users can belong to one or more groups, which determine
whether they can peform certain actions; adding a user to the docker group allows them to execute docker commands).  To add the user with
usrname `<username>`, to the docker group, run:
```
usermod -aG docker <username>
```

> _Some Linux systems may not have_ `usermod` _installed, as there are different programs that can be used to change user settings;_
> _please consult your Linux distribution documentation if necessary._

#### Installation

There are a number of ways to install Companion; details for an installation using Docker are described below. If you encounter an issue when installing Companion please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/companion/issues) or email us at path-help@sanger.ac.uk.

The easiest way to use the pipeline is to use the prepared [Docker image](https://hub.docker.com/r/sangerpathogens/companion/) which contains all external dependencies.
```
docker pull sangerpathogens/companion
```

#### Usage

##### Local copy of Companion

To create a local copy of companion, you can download this repo from github (if you are familiar with github, you may
of course prefer to _clone_ or _fork_ it).
```
curl -L -o companion-master.zip https://github.com/sanger-pathogens/companion/archive/master.zip  # or click the green button on the guthub web page
unzip companion-master.zip
mv companion-master my-companion-project # renaming it to something meaningful to you is a good idea
```

Now you can run Companion.   There is an example dataset and parameterization included in the distribution, so
to get started just run:
```
nextflow run my-companion-project -profile docker
```
The argument `-profile docker` instructs nextflow to run the sangerpathogens/companion docker image for the dependencies;
the `nextflow.config` file (and files referenced within it) define the docker profile and the docker image to be used.

##### Running Companion direct from a repository

If you run nextflow with the name of a github repository, it will pull (download) the contents of the repository and run with those.
For example, the following command will do the same as the "local copy" example above:
```
nextflow run sanger-pathogens/companion -profile docker
```
It is best to use this with some caution.  After the command above is
run, nextflow will have stored a local copy of the repository in `.nextflow/assets/sanger-pathogens`
(note that `.nextflow` is a hidden directory, and will not usually be visible; use the command `ls -la .nextflow` to see it).

If you run the same command again it will this time use the _local_ copy instead of pulling a copy from the repository.  You can
edit the files in your local copy, and nextflow will work from your (now different) version of sanger-pathogens/companion.

If you are familiar with repositories, and the workflow appropriate to using them, this can be a very convenient way of
working.  You can create your own github repository to store and share your work, and track versions.

If you are not familiar with git repositories, it can become quite confusing, and you should probably work with a simple local copy.


##### Preparing reference annotations

Further documentation on preparing reference data can be found in the
[GitHub wiki](https://github.com/sanger-pathogens/companion/wiki/Preparing-reference-data-sets).


## License
Companion is free software, licensed under [ISC](https://github.com/sanger-pathogens/companion/blob/master/LICENSE).

## Feedback/Issues
Please report any issues to the [issues page](https://github.com/sanger-pathogens/companion/issues) or email path-help@sanger.ac.uk

## Citation
If you use this software please cite:
__Companion: a web server for annotation and analysis of parasite genomes.__
Steinbiss S, Silva-Franco F, Brunk B, Foth B, Hertz-Fowler C et al.
Nucleic Acids Research, 44:W29-W34, 2016.  
DOI: [10.1093/nar/gkw292](http://dx.doi.org/10.1093/nar/gkw292)
