[![Snakemake](https://img.shields.io/badge/snakemake-≥5.5.0-brightgreen.svg?style=flat-square)](https://snakemake.bitbucket.io)
[![License](https://img.shields.io/apm/l/vim-mode.svg)](LICENSE)

# varCA
A pipeline for running an ensemble of variant callers to predict variants in ATAC-seq reads.

The entire pipeline is made up of two smaller pipelines. The `prepare` pipeline calls each variant caller and prepares the resulting data for use by the `classify` pipeline, which runs the ensemble classifier to predict the existence of variants at each site.

# execution
```
# install snakemake via conda (if not already installed)
conda install -c bioconda -c conda-forge -n snakemake 'snakemake>=5.5.0'
# activate the conda env in which it was installed
conda activate snakemake
# execute the pipeline on example data on an SGE cluster
qsub run.bash
```

The pipeline is written as Snakefiles, so it must be executed via [Snakemake](https://snakemake.readthedocs.io/en/stable/). See the [`run.bash` script](run.bash) for an example. Make sure to provide required input and options in the [config files](configs) before executing.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a BAM). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipelines.

# dependencies
We highly recommend you install [Snakemake via conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda) so that you can use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipelines. Otherwise, you must manually install the dependencies listed in the [env files](envs).

# files and directories

### [Snakefiles/Snakefile-prepare](Snakefiles/Snakefile-prepare)
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for preparing data for the classifier. It uses ATAC-seq FASTQ files to generate a tab-delimited table containing variant caller output for every site in open chromatin regions of the genome. The output of this pipeline can be fed into `Snakefiles/Snakefile-classify`.

### [Snakefiles/Snakefile-classify](Snakefiles/Snakefile-classify)
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for training and testing the classifier. It uses the output of `Snakefiles/Snakefile-prepare`.

### [configs/](configs)
Config files that define options and input for the `prepare` and `classify` pipelines. You should start by filling these out.

### [callers/](callers)
Scripts for executing each of the variant callers which are used by the `prepare` pipeline. Small pipelines can be written for each caller by using a special naming convention. See the [caller README](callers/README.md) for more information.

### [breakCA/](breakCA)
Scripts for calculating posterior probabilities for the existence of an insertion or deletion, which can be used when running the classifier. These scripts are an adaptation from [@Arkosen](https://github.com/Arkosen)'s [breakCA code](https://github.com/Arkosen/BreakCA/tree/master/bin).

### [scripts/](scripts)
Various scripts used by the pipeline. See the [script README](scripts/README.md) for more information.

### [run.bash](run.bash)
An example bash script for executing the `prepare` and `classify` pipelines on an SGE cluster using `snakemake` and `conda`.
