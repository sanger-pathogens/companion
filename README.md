# Companion
A portable, scalable eukaryotic genome annotation pipeline implemented in Nextflow.

[![Build Status](https://travis-ci.org/sanger-pathogens/companion.svg?branch=master)](https://travis-ci.org/sanger-pathogens/companion)  
[![License: ISC](https://img.shields.io/badge/License-ISC-brightgreen.svg)](https://github.com/sanger-pathogens/companion/blob/master/LICENSE)  
[![status](https://img.shields.io/badge/NAR-10.1093%2Fnar.gkw292-brightgreen.svg)](https://doi.org/10.1093/nar/gkw292)

## Content
  * [Introduction](#introduction)
  * [Installation](#installation)
    * [Required dependencies](#required-dependencies)
  * [Usage](#usage)
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

## Installation
Companion has the following dependencies:

### Required dependencies
 * [Nextflow](http://nextflow.io)

To install Nextflow, run:
```
curl -fsSL get.nextflow.io | bash
```

There are a number of ways to install Companion. Details for installing using Docker are described below. If you encounter an issue when installing Companion please contact your local system administrator. If you encounter a bug please log it [here](https://github.com/sanger-pathogens/companion/issues) or email us at path-help@sanger.ac.uk.

The easiest way to use the pipeline is to use the prepared [Docker container](https://hub.docker.com/r/sangerpathogens/companion/) which contains all external dependencies.
```
docker pull sangerpathogens/companion
```
## Usage

Start an example run using Docker (using the example dataset and parameterization included in the distribution):
```
nextflow run sanger-pathogens/companion -profile docker
```

For your own runs, provide your own file names, paths, parameters, etc. as defined in the `nextflow.config` file.

### Preparing reference annotations

The reference annotations used in the pipeline need to be pre-processed before they can be used. See the the [GitHub wiki](https://github.com/sanger-pathogens/companion/wiki/Preparing-reference-data-sets) for more details. There are also pre-generated reference sets for various parasite species/families.

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