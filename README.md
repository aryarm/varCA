[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.0-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# breakCA
A pipeline for running an ensemble of variant callers to predict variants in ATAC-seq reads.

This project is based on [@Arkosen](https://github.com/Arkosen)'s [project of the same name](https://github.com/Arkosen/BreakCA).

# files and directories

### Snakefiles/Snakefile-prepare
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for preparing data for the classifier. It uses ATAC-seq FASTQ files to generate a tab-delimited table containing variant caller output for every site in open chromatin regions of the genome. The output of this pipeline can be fed into `Snakefiles/Snakefile-classify`.

### Snakefiles/Snakefile-classify
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for training and testing the classifier. It uses the output of `Snakefiles/Snakefile-prepare`.

### configs/
Config files that define options and input for the `prepare` and `classify` pipelines.

### callers/
Scripts for executing each of the variant callers. Small pipelines can be written for each caller by using a special naming convention. See the [caller README](https://github.com/aryam7/breakCA/blob/master/callers/README.md) for more information.

### breakCA/
Scripts for calculating posterior probabilities for the existence of an insertion or deletion, which can be used when prearing data for the classifier. These scripts are an adaptation from [@Arkosen](https://github.com/Arkosen)'s [breakCA code](https://github.com/Arkosen/BreakCA/tree/master/bin).

### scripts/
Various scripts used by the pipeline. See the [script README](https://github.com/aryam7/breakCA/blob/master/callers/README.md) for more information.

### run-prepare.bash
An example bash script for executing the `prepare` pipeline on an SGE cluster using snakemake.

### run-classify.bash
An example bash script for executing the `classify` pipeline on an SGE cluster using snakemake.

# execution
The pipeline is written as Snakefiles and so can be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/). See the `run-<pipeline>.bash` scripts for an example. Make sure to provide required input and options in the `config files` before executing.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a BAM). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipeline.

# dependencies
If you have [conda](https://conda.io/docs/user-guide/install/index.html) installed (highly recommended), use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipelines. Otherwise, you must manually install the dependencies listed in the [env files](https://github.com/aryam7/breakCA/blob/master/envs).
