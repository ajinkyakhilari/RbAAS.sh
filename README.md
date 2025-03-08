# RbAAS.sh
Reference based Amplicon Assembly Script

# Multi-Environment Variant Calling Pipeline

This repository contains a Bash pipeline for processing sequencing samples—merging FASTQ files, filtering reads, aligning to a reference genome, trimming primers, performing variant calling with Clair3, and generating a consensus sequence with masking. Each software tool used in the pipeline is installed in its own isolated Conda environment defined by separate YAML files.

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [Installation](#installation)
  - [Conda Installation](#conda-installation)
  - [Creating Environments](#creating-environments)
- [Pipeline Usage](#pipeline-usage)
- [Automatic Environment Creation](#automatic-environment-creation)
- [Environment YAML Files](#environment-yaml-files)
- [License](#license)

## Overview

This pipeline is designed for Oxford Nanopore sequencing data. It includes the following steps:

1. **Quality Filtering:** Uses [fastp](https://github.com/OpenGene/fastp) to filter reads.
2. **Alignment:** Uses [minimap2](https://github.com/lh3/minimap2) and [samtools](http://www.htslib.org/) for alignment and sorting.
3. **Primer Trimming:** Uses [ivar](https://github.com/andersen-lab/ivar) to trim primer sequences.
4. **Variant Calling:** Uses [Clair3](https://github.com/HKU-BAL/Clair3) for variant calling.
5. **Masking and Consensus:** Uses [bcftools](http://www.htslib.org/doc/bcftools.html) and [tabix](https://www.htslib.org/doc/tabix.html) to mask low-confidence regions and generate a final consensus.

## Requirements

- [Conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda or Anaconda)
- Git

## Installation

### Conda Installation

If you don't have Conda installed, please follow the instructions for [Miniconda](https://docs.conda.io/en/latest/miniconda.html).

### Creating Environments

This repository provides separate YAML files for each required tool:

- **fastp:** `fastp_env.yml`
- **minimap2:** `minimap2_env.yml`
- **samtools:** `samtools_env.yml`
- **ivar:** `ivar_env.yml`
- **bcftools:** `bcftools_env.yml`
- **tabix (htslib):** `tabix_env.yml`
- **Clair3:** `clair3_env.yml`

To create each environment, run:

```bash
conda env create -f fastp_env.yml
conda env create -f minimap2_env.yml
conda env create -f samtools_env.yml
conda env create -f ivar_env.yml
conda env create -f bcftools_env.yml
conda env create -f tabix_env.yml
conda env create -f clair3_env.yml

