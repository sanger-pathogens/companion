# Companion
A portable, scalable eukaryotic genome annotation pipeline implemented in Nextflow.

[![Build Status](https://travis-ci.org/sanger-pathogens/companion.svg?branch=master)](https://travis-ci.org/sanger-pathogens/companion)  
[![License: ISC](https://img.shields.io/badge/License-ISC-brightgreen.svg)](https://github.com/sanger-pathogens/companion/blob/master/LICENSE)  
[![status](https://img.shields.io/badge/NAR-10.1093%2Fnar.gkw292-brightgreen.svg)](https://doi.org/10.1093/nar/gkw292)


## Content
 * [Introduction](#introduction)
 * [Dependencies](#dependencies)
   * [Docker](#docker)
 * [Installation](#installation)
 * [Usage](#usage)
   * [Local copy of Companion](#local-copy-of-companion)
   * [Running Companion direct from a repository](#running-companion-direct-from-a-repository)
   * [Preparing reference annotations](#preparing-reference-annotations)
 * [License](#license)
 * [Feedback/Issues](#feedbackissues)
 * [Citation](#citation)


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

## Dependencies

Companion has the following dependencies:

 * Java 8 or later
 * [Nextflow](http://nextflow.io)
 * [Docker](https://www.docker.com/) (if using the Docker image to satisfy dependencies)

To check if you have Java installed, and the version, use the command `java -version`.  Note that this will give you a version number
of 1.8 for Java 8, 1.9 for Java 9, etc.

If you need to install Java, on an Ubuntu system run:
```
apt-get install default-jre
```
> _For other Linux systems, please consult your distribution documentation._

To install Nextflow, run:
```
curl -fsSL get.nextflow.io | bash
```
This will create an executable called 'nextflow', which should be moved to a suitable directory, for example:
```
mv nextflow /usr/local/bin/
```
Use the command `which nextflow` to check that it is found in your path.

### Docker

Docker is required if you intended to use the Docker image, as recommended below, to satisfy the dependencies.

To install Docker, see the installation guide for
[Ubuntu](https://docs.docker.com/install/linux/docker-ce/ubuntu/),
[Centos](https://docs.docker.com/install/linux/docker-ce/centos/),
[Debian](https://docs.docker.com/install/linux/docker-ce/debian/) or
[Fedora](https://docs.docker.com/install/linux/docker-ce/fedora/).

Users running Companion with Docker will need to be added to the `docker` group (unix users can belong to one or more groups, which determine
whether they can peform certain actions; adding a user to the docker group allows them to execute docker commands).  To add user `<username>`, to
the docker group, run:
```
usermod -aG docker <username>
```
> _Some Linux systems may not have_ `usermod` _installed, as there are different programs that can be used to change user settings;_
> _please consult your Linux distribution documentation if necessary._

## Installation

There are a number of ways to install Companion. Details for an installation using Docker are described below. If you encounter an issue when installing Companion please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/companion/issues) or email us at path-help@sanger.ac.uk.

The easiest way to use the pipeline is to use the prepared [Docker image](https://hub.docker.com/r/sangerpathogens/companion/) which contains all external dependencies.
```
docker pull sangerpathogens/companion
```

## Usage

### Local copy of Companion

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
The argument `-profile docker` instructs nextflow to run the sangerpathogens/companion docker image for the dependencies.

Have a look at the `nextflow.config` file to see the definition of the docker profile, and how the docker image is specified.
You will also find file names, paths, parameters, etc. that you can edit to perform your own runs.  The following warrant
a special mention:

*inseq*  The input FASTA file  (`${baseDir}/example-data/L_donovani.1.fasta` in the example parameter file included wirth the distribution)

*ref_dir* The directory containing reference genomes (`${baseDir}/example-data/references` in the example file)

*dist_dir* The directory that will contain the newly created output files (`${baseDir}/example-data-output` in the example file)

*run_snap* We recommend SNAP is disabled, as it has not provided useful results in this pipeline (`false` in the example file)


### Running Companion direct from a repository

If you run nextflow with the name of a github repository, it will pull the contents of the repository and run with those.
This command will do the same as the "local copy" example above:
```
nextflow run sanger-pathogens/companion -profile docker
```
It is best to use this with some caution.  After the command above is
run, nextflow will have stored a local copy of the repository in `.nextflow/assets/sanger-pathogens`, and if you run
the command again it will this time use the _local_ copy instead of pulling a copy from the repository.  You can
edit the files in your local copy, and nextflow will work from your (now different) version of sanger-pathogens/companion.

If you are familiar with repositories, and the workflow appropriate to using them, this can be a very convenient way of
working;   otherwise it can become quite confusing, and you may find it easier to work with a simple local copy.

### Preparing reference annotations

The reference annotations used in the pipeline need to be pre-processed before they can be used.  Only a few pre-generated
reference sets for various parasite species/families are included in the distribution as examples.

To add a reference organism, you will need:

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
1. Add new section to `amber-test-data/references/references-in.json`, using the
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

You can now run _Companion_, and the new reference will be included.

Further documentation on preparing reference data can be found in the [GitHub wiki](https://github.com/sanger-pathogens/companion/wiki/Preparing-reference-data-sets).


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
