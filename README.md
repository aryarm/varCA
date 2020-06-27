[![Snakemake](https://img.shields.io/badge/snakemake-5.18.0-brightgreen.svg?style=flat-square)](https://snakemake.readthedocs.io/)
[![License](https://img.shields.io/apm/l/vim-mode.svg)](LICENSE)

# varCA
A pipeline for running an ensemble of variant callers to predict variants in ATAC-seq reads.

The entire pipeline is made up of two smaller subworkflows. The `prepare` subworkflow calls each variant caller and prepares the resulting data for use by the `classify` subworkflow, which runs the ensemble classifier to predict the existence of variants at each site.

# download
Execute the following commands or download the [latest release](https://github.com/aryam7/varCA/releases/latest) manually.
```
wget -O- -q https://github.com/aryam7/varca/tarball/master | tar mxvzf -
mv aryam7-* varca
```
Also consider downloading the [example data](https://github.com/aryam7/varCA/releases/latest/download/data.tar.gz).
```
cd varca
wget -O- -q https://github.com/aryam7/varCA/releases/latest/download/data.tar.gz | tar xvzf -
```
The example data includes the following files:

- `samples.tsv` - An example samples file; required for the master pipeline and the `predict` subworkflow
- `snp.rda` - A trained RF model for classifying SNVs
- `indel.rda` - A trained RF model for classifying indels
- `snp.tsv.gz` - Prepared, example SNV training data
- `indel.tsv.gz` - Prepared, example indel training data
- `hg38.chr1.*` - Chromosome 1 of the hg38 reference genome and its index files
- `jurkat.chr1.bam` - BAM file from ATAC-seq on chromosome 1 from the Jurkat sample in GSE129086
- `jurkat.chr1.bam.bai` - Index for the Jurkat BAM file
- `jurkat.chr1.bed` - ATAC-seq peaks for the Jurkat sample
- `molt4.chr1.*.fq.gz` - FASTQ files from paired ATAC-seq reads for chromosome 1 from the Molt-4 sample in GSE129086

# execution
On example data:
```
# 1) install snakemake via conda (if not already installed)
conda install -c bioconda -c conda-forge 'snakemake==5.18.0'

# 2) execute the pipeline on example data locally
./run.bash &

# OR execute the pipeline on example data on an SGE cluster
qsub run.bash
```

The pipeline is written as a Snakefile, so it must be executed via [Snakemake](https://snakemake.readthedocs.io). See the [`run.bash` script](run.bash) for an example. Make sure to provide required input and options in the [config files](configs) before executing. The config files are currently configured to run the pipeline on the example data provided.

### If this is your first time using Snakemake
We highly recommend that you run `snakemake --help` to learn about all of the options available to you. You might discover, for example, that calling Snakemake with the `-n -p -r` flags can be a helpful way to check that the pipeline will be executed correctly before you run it. This can also be a good way to familiarize yourself with the steps of the pipeline and their inputs and outputs (the latter of which are inputs to the first rule in each workflow -- ie the `all` rule).

Another important thing to know is that Snakemake will not recreate output that it has already generated, unless you request it. If a job fails or is interrupted, subsequent executions of Snakemake will just pick up where it left off. This can also apply to files that *you* create and provide in place of the files it would have generated.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a BAM). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipeline.

# dependencies
We highly recommend you install [Snakemake via conda](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda) so that you can use the `--use-conda` flag when calling `snakemake` to let it automatically handle all dependencies of the pipeline. Otherwise, you must manually install the dependencies listed in the [env files](envs).
Although currently untested, we also provide the option of executing the pipeline in a Docker container using [`singularity`](https://sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps). Just provide Snakemake with the `--use-conda --use-singularity` flags when you execute it. Unlike with the previous method, having `conda` installed on your machine is not a requirement for this option, since it will be installed automatically in the Docker container.

# files and directories

### [Snakefile](Snakefile)
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for calling variants from a set of ATAC-seq reads. This pipeline is made up of two subworkflows:

1. the [`prepare` subworkflow](rules/prepare.smk), which prepares the reads for classification and
2. the [`classify` subworkflow](rules/classify.smk), which creates a VCF containing predicted variants

### [rules/](rules)
Snakemake rules for the `prepare` and `classify` subworkflows. You can either execute these subworkflows from the [master Snakefile](#snakefile) or individually as their own Snakefiles. See the [rules README](rules/README.md) for more information.

### [configs/](configs)
Config files that define options and input for the pipeline and the `prepare` and `classify` subworkflows. If you want to predict variants from your own ATAC-seq data, you should start by filling out [the config file for the pipeline](/configs#configyaml).

### [callers/](callers)
Scripts for executing each of the variant callers which are used by the `prepare` subworkflow. Small pipelines can be written for each caller by using a special naming convention. See the [caller README](callers/README.md) for more information.

### [breakCA/](breakCA)
Scripts for calculating posterior probabilities for the existence of an insertion or deletion, which can be used as features for the classifier. These scripts are an adaptation from [@Arkosen](https://github.com/Arkosen)'s [BreakCA code](https://www.biorxiv.org/content/10.1101/605642v1.abstract).

### [scripts/](scripts)
Various scripts used by the pipeline. See the [script README](scripts/README.md) for more information.

### [run.bash](run.bash)
An example bash script for executing the pipeline using `snakemake` and `conda`. Any parameters to this script are passed directly to `snakemake`.
