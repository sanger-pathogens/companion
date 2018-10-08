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

#### Check

. `java -version` should say you have Java 1.8 or greater
. `nextflow info` will print system information if nextflow has been installed successfully
. `systemctl status docker` will tell you if docker is active (running)

### 2. Install _Companion_

Execute these command in the directory you want to keep your _Companion_ work in.  Do this as a normal user, i.e. not as root or using sudo.
```
curl -L -o companion-master.zip https://github.com/sanger-pathogens/companion/archive/master.zip
unzip companion-master.zip
mv companion-master my-companion-project # renaming it to something meaningful to you is a good idea
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

*ABACAS_BIN_CHR* <add description> (`LDON_0` in the example file)

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

<adapt from README>

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
