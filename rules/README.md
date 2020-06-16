## The `prepare` subworkflow
The [`prepare` subworkflow](prepare.smk) is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for preparing data for the classifier. It uses ATAC-seq FASTQ (or BAM) files to generate a tab-delimited table containing variant caller output for every site in open chromatin regions of the genome. The `prepare` subworkflow uses the scripts in the [callers directory](callers) to run every variant caller in the ensemble.

The `prepare` subworkflow is included within the [master pipeline](/Snakefile) automatically. However, you can also execute the `prepare` subworkflow on its own, as a separate Snakefile.
First, make sure that you fill out the [`prepare.yaml`](/configs/prepare.yaml) and the [`callers.yaml`](/configs/callers.yaml) config files, which specify input and parameters for the `prepare` subworkflow. See the [config README](/configs) for more information.
You can execute the workflow like this:

	snakemake -s rules/prepare.smk --use-conda -j

The output of the `prepare` pipeline will be in `<output_directory>/merged_<variant_type>/<sample_ID>/final.tsv.gz`

## The `classify` subworkflow
The [`classify` subworkflow](classify.smk) is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for training and testing the classifier. It uses the TSV output by `prepare` subworkflow. Its final output is a VCF containing predicted variants.

The `classify` subworkflow is included within the [master pipeline](/Snakefile) automatically. However, you can also execute the `classify` subworkflow on its own, as a separate Snakefile.
First, make sure that you fill out the [`classify.yaml`](/configs/classify.yaml) config file, which specifies input and parameters for the `classify` subworkflow. See the [config README](/configs) for more information.
You can execute the workflow like this:

	snakemake -s rules/classify.smk --use-conda -j

The output of the `classify` pipeline will be in `<output_directory>/classify/<sample_ID>/final.vcf.gz`

## Creating your own trained model
You may want to create your own trained models (rather than use the ones we provided in the example data) for any number of reasons. The most common are

1. You changed a parameter in the `callers.yaml` config file, or
2. You'd like to include a new variant caller that hasn't already been implemented as a caller script in the [callers directory](/callers)

For the sake of this example, let's say you'd like to include a new indel variant caller (ie #2). You've also already followed the directions in the [callers README](/callers/README.md) to create your own caller script, and you've modified the `prepare.yaml` and `callers.yaml` config files to include your new indel caller. However, before you can predict variants using the indel caller, you must create a new trained classification model that knows how to interpret this new input.

To do this, we recommend first downloading the training and truth sets we used to create our model. First, download the [GM12878 FASTQ files from Buenrostro et al](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753). Then, download the corresponding [Platinum Genomes VCF for that sample](https://www.illumina.com/platinumgenomes.html).

Fill out the `prepare.yaml` config file with the necessary information to run it on the GM12878 sample. Add `pg-indel` as the last caller under the `indel_callers` list. Also, include the path to your Platinum Genomes VCF in the `callers.yaml` config file under `pg`, `pg-snp`, and `pg-indel`. Next, execute the `prepare` subworkflow on its own.

After the `prepare` subworkflow has finished running, add the sample (add specfically, the path to its `final.tsv.gz` file) as a dataset in the `classify.yaml` file. Most of the yaml for this has already been written for you. Next, execute the `classify` subworkflow on its own.

We recommend testing your trained model before using it. You can include test sets as datasets in the `classify.yaml` file. For example, you could test on odd chrosomes by subsetting the `final.tsv.gz` file and providing that as a separate "predict" dataset.
