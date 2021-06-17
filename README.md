[![Snakemake](https://img.shields.io/badge/snakemake-5.18.0-brightgreen.svg?style=flat-square)](https://snakemake.readthedocs.io/)
[![License](https://img.shields.io/apm/l/vim-mode.svg)](LICENSE)

# varCA
A pipeline for running an ensemble of variant callers to predict variants from ATAC-seq reads.

The entire pipeline is made up of two smaller subworkflows. The `prepare` subworkflow calls each variant caller and prepares the resulting data for use by the `classify` subworkflow, which uses an ensemble classifier to predict the existence of variants at each site.

### [Code Ocean](https://codeocean.com/capsule/6980349/tree/v1)
Using [our Code Ocean compute capsule](https://codeocean.com/capsule/6980349/tree/v1), you can execute [VarCA v0.2.1](https://github.com/aryarm/varCA/releases/tag/v0.2.1) on example data without downloading or setting up the project. To interpret the output of VarCA, see the output sections of the [`prepare` subworkflow](rules#output) and the [`classify` subworkflow](rules#output-1) in the [rules README](rules/README.md).

# download
Execute the following command or download the [latest release](https://github.com/aryarm/varCA/releases/latest) manually.
```
git clone https://github.com/aryarm/varCA.git
```
Also consider downloading the [example data](https://github.com/aryarm/varCA/releases/latest/download/data.tar.gz).
```
cd varCA
wget -O- -q https://github.com/aryarm/varCA/releases/latest/download/data.tar.gz | tar xvzf -
```

# setup
The pipeline is written as a Snakefile which can be executed via [Snakemake](https://snakemake.readthedocs.io). For most users, we recommend installing Snakemake via mamba as [described in their documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html#installation-via-conda-mamba).

However, if your aim is complete reproducbility, we recommend installing version 5.24.2 via this command:
```
conda create -n snakemake -c bioconda -c conda-forge --no-channel-priority 'snakemake==5.24.2'
```

# execution
1. Activate snakemake via `conda`:
    ```
    conda activate snakemake
    ```
2. Execute the pipeline on the example data

    Locally:
    ```
    ./run.bash &
    ```
    __or__ on an SGE cluster:
    ```
    ./run.bash --sge-cluster &
    ```
#### Output
VarCA will place all of its output in a new directory (`out/`, by default). Log files describing the progress of the pipeline will also be created there: the `log` file contains a basic description of the progress of each step, while the `qlog` file is more detailed and will contain any errors or warnings. You can read more about the pipeline's output in the [rules README](rules/README.md).

#### Executing the pipeline on your own data
You must modify [the config.yaml file](configs#configyaml) to specify paths to your data. The config file is currently configured to run the pipeline on the example data provided.

#### Executing each portion of the pipeline separately
The pipeline is made up of [two subworkflows](rules). These are usually executed together automatically by the master pipeline, but they can also be executed on their own for more advanced usage. See the [rules README](rules/README.md) for execution instructions and a description of the outputs. You will need to execute the subworkflows separately [if you ever want to create your own trained models](rules#training-and-testing-varca).

#### Reproducing our results
We provide the example data so that you may quickly (in ~1 hr, excluding dependency installation) verify that the pipeline can be executed on your machine. This process does not reproduce our results. Those with more time can follow [these steps](rules#testing-your-model--reproducing-our-results) to create all of the plots and tables in our paper.

### If this is your first time using Snakemake
We recommend that you run `snakemake --help` to learn about Snakemake's options. For example, to check that the pipeline will be executed correctly before you run it, you can call Snakemake with the `-n -p -r` flags. This is also a good way to familiarize yourself with the steps of the pipeline and their inputs and outputs (the latter of which are inputs to the first rule in each workflow -- ie the `all` rule).

Note that Snakemake will not recreate output that it has already generated, unless you request it. If a job fails or is interrupted, subsequent executions of Snakemake will just pick up where it left off. This can also apply to files that *you* create and provide in place of the files it would have generated.

By default, the pipeline will automatically delete some files it deems unnecessary (ex: unsorted copies of a BAM). You can opt to keep these files instead by providing the `--notemp` flag to Snakemake when executing the pipeline.

# files and directories

### [Snakefile](Snakefile)
A [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for calling variants from a set of ATAC-seq reads. This pipeline automatically executes two subworkflows:

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
An example bash script for executing the pipeline using `snakemake` and `conda`. Any arguments to this script are passed directly to `snakemake`.
