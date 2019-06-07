[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.0-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)

# merge_callers
A Snakemake pipeline for running multiple variant callers and merging their output into a single, large table.

# files

### Snakefile
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline defining rules for every step of the analysis. It uses DNA and RNA FASTQ files to generate a summary of allelic imbalance for each gene.

### config.yaml
Defines options and input for the Snakemake pipeline.

### run-all.bash
An example bash script for executing the entire pipeline on an SGE cluster using snakemake.

# execution
The pipeline is written as Snakefiles and so can be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/). See the [run-all.bash](https://github.com/aryam7/as_analysis/blob/master/run-all.bash) script for an example. Make sure to provide required input and options in the [config file](https://github.com/aryam7/as_analysis/blob/master/config.yaml) before executing.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a BAM). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipeline.

# dependencies
If you have [conda](https://conda.io/docs/user-guide/install/index.html) installed (highly recommended), use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipeline. Otherwise, you must manually install the dependencies listed in the [env.yml](https://github.com/aryam7/merge_callers/blob/master/env.yml) file.
