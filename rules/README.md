# The varCA pipeline
![Pipeline Skeleton](https://drive.google.com/uc?export=view&id=10a8wgPMPojUR5o15_Ik9TOsAvVg63HK2)

The pipeline consists of the `prepare` (red) and `classify` (green) subworkflows. The prepare subworkflow runs a set of variant callers on the provided ATAC-seq data and prepares their output for use by the `classify` subworkflow, which uses a trained ensemble classifier to predict the existence of variants.
The `prepare` subworkflow can use FASTQ or BAM/BED files as input. The `classify` subworkflow can predict variants from the prepared datasets or use them for training/testing.
If a pre-trained model is available (orange), the two subworkflows can be executed together automatically via the master pipeline. However the subworkflows must be executed separately for training and testing (see [below](#training-and-testing-varca)).

## The `prepare` subworkflow
The [`prepare` subworkflow](prepare.smk) is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for preparing data for the classifier. It generates a tab-delimited table containing variant caller output for every site in open chromatin regions of the genome. The `prepare` subworkflow uses the scripts in the [callers directory](callers) to run every variant caller in the ensemble.

### execution
The `prepare` subworkflow is included within the [master pipeline](/Snakefile) automatically. However, you can also execute the `prepare` subworkflow on its own, as a separate Snakefile.

First, make sure that you fill out the [`prepare.yaml`](/configs/prepare.yaml) and the [`callers.yaml`](/configs/callers.yaml) config files, which specify input and parameters for the `prepare` subworkflow. See the [config README](/configs) for more information.

Then, just call Snakemake with `-s rules/prepare.smk`:

	snakemake -s rules/prepare.smk --use-conda -j

### output
The primary outputs of the `prepare` pipeline will be in `<output_directory>/merged_<variant_type>/<sample_ID>/final.tsv.gz`. However, several intermediary directories and files are also generated:

- `align/` - output from the BWA FASTQ alignment step and samtools PCR duplicate removal steps
- `peaks/` - output from MACS 2 and other files required for calling peaks
- `callers/` - output from each [caller script](/callers) in the ensemble (see the [callers README](/callers/README.md) for more information) and the variant normalization and feature extraction steps
- `merged_<variant_type>/` - all other output in the `prepare` subworkflow, including the merged and final datasets for each variant type (ie SNV or indels)

## The `classify` subworkflow
The [`classify` subworkflow](classify.smk) is a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for training and testing the classifier. It uses the TSV output from the `prepare` subworkflow. Its final output is a VCF containing predicted variants.

### execution
The `classify` subworkflow is included within the [master pipeline](/Snakefile) automatically. However, you can also execute the `classify` subworkflow on its own, as a separate Snakefile.

First, make sure that you fill out the [`classify.yaml`](/configs/classify.yaml) config file, which specifies input and parameters for the `classify` subworkflow. See the [config README](/configs) for more information.

Then, just call Snakemake with `-s rules/classify.smk`:

	snakemake -s rules/classify.smk --use-conda -j

### output
The primary output of the `classify` subworkflow will be in `<output_directory>/classify/<sample_ID>/final.vcf.gz`. However, several intermediary files are also generated:

#### for all datasets

- `subset.tsv.gz` - the result of subsetting some callers from the input (only if requested)
- `filter.tsv.gz` - the result of filtering rows from the input (only if requested)
- `prepared.tsv.gz` - datasets that are ready to be interpreted by the ensemble classifier

#### for prediction datasets

- `results.tsv.gz` - the predicted variants in TSV format
- `results.vcf.gz` - the predicted variants in VCF format with recalibrated QUAL scores

#### for training datasets

- `model.rda` - the trained classifier
- `variable_importance.tsv` - a table containing importance values for every feature in the RF classifier; visualize this with [`importance_plot.py`](/scripts/importance_plot.py)
- `tune_matrix.tsv` - the results from hyperparameter tuning (only if requested)
- `tune.pdf` - a visualization of the results in tune_matrix.tsv

#### for test datasets

- `prc/results.pdf` - precision-recall plot for evaluating the performance of varCA and comparing it to other variant callers
- `prc/curves` - data used to create the precision-recall plot
- `prc/pts` - single point metrics evaluating the performance of varCA and comparing it to other callers; you may aggregate these with [`metrics_table.py`](/scripts/metrics_table.py)

# Training and Testing varCA
## Motivation
### Training
You may want to create your own trained models (rather than use the ones we provided in the example data) for any number of reasons. The most common are

1. You changed one of the caller specific parameters in the `callers.yaml` config file, or
2. You would like to use a different set of variant callers than the example models support, or
3. You'd like to include a new variant caller that hasn't already been implemented as a caller script in the [callers directory](/callers)
4. You'd like to wholly and completely reproduce our results (since the training and testing steps are usually skipped by the master pipeline)

In general, you will need to create a new trained model whenever the columns in the output of the `prepare` subworkflow change.
### Testing
You may want to test a trained model if:

1. You want to create precision-recall curves and generate performance metrics for every variant caller in the ensemble
2. You want to compare the performance of varCA to other variant callers
3. You'd like to reproduce our results (since the training and testing steps are usually skipped by the master pipeline)

## Creating your own trained model
For the sake of this example, let's say you'd like to include a new indel variant caller (ie #3 above). You've also already followed the directions in the [callers README](/callers/README.md) to create your own caller script, and you've modified the `prepare.yaml` and `callers.yaml` config files to include your new indel caller. However, before you can predict variants using the indel caller, you must create a new trained classification model that knows how to interpret your new input.

To do this, we recommend downloading the truth set we used to create our model. First, download the [GM12878 FASTQ files from Buenrostro et al](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753). Specify the path to these files in `data/samples.tsv`, the samples file within the example data. Then, download the corresponding [Platinum Genomes VCF for that sample](https://www.illumina.com/platinumgenomes.html).

Fill out the `prepare.yaml` config file with the necessary information to run it on the GM12878 sample. Add `pg-indel` as the last caller under the `indel_callers` list, so that the `prepare` subworkflow will know to include the Platinum Genomes VCF in its output. Also, include the path to your Platinum Genomes VCF in the `callers.yaml` config file under `pg`, `pg-snp`, and `pg-indel`. Next, execute the `prepare` subworkflow on its own (see [above](#execution)).

After the `prepare` subworkflow has finished running, add the sample (specfically, the path to its `final.tsv.gz` file) as a dataset in the `classify.yaml` file and specify its attributes. Most of the yaml for this has already been written for you. Specifically, make sure you've specified `pg-indel` (for the Platinum Genomes VCF) as the truth caller ID. Next, execute the `classify` subworkflow on its own (see [above](#execution-1)).

The `classify` subworkflow can only create one trained model at a time, so you will need to repeat these steps if you'd also like to create a trained model for SNVs. Just replace every mention of "indel" in `classify.yaml` with "snp". Also remember to use only the SNV callers (ie GATK, VarScan 2, and VarDict).

## Testing your model / Reproducing our Results
For this example, we will demonstrate how you can reproduce the results in our paper using the `indel.tsv.gz` truth dataset we provided in the example data. This data was generated by running the `prepare` subworkflow on the GM12878 data as described [above](#creating-your-own-trained-model). If you [ran the `prepare` subworkflow to create your own trained model](#creating-your-own-trained-model), just use your truth dataset instead of the one we provided in the example data.

First, split the truth dataset by chromosome parity using `awk` commands like this:

	zcat data/indel.tsv.gz | { read -r head && echo "$head" && awk -F $"\t" -v 'OFS=\t' '$1 ~ /^[0-9]+$/ && $1%2'; } | gzip > data/indel-odds.tsv.gz
	zcat data/indel.tsv.gz | { read -r head && echo "$head" && awk -F $"\t" -v 'OFS=\t' '$1 ~ /^[0-9]+$/ && !($1%2)'; } | gzip > data/indel-evens.tsv.gz

Then, specify in `classify.yaml` the path to the `indel-odds.tsv.gz` file as your training set and the path to the `indel-evens.tsv.gz` file as your test set. Make sure that both your training and test sets have a truth caller ID and that you've specified the test set ID in the list of `predict` datasets. Then, execute the `classify` subworkflow on its own (see [above](#execution-1)).

The `classify` subworkflow can only work with one trained model at a time, so you will need to repeat these steps if you'd also like to create a trained model for SNVs. Just replace every mention of "indel" in `classify.yaml` with "snp". Also remember to use only the SNV callers (ie GATK, VarScan 2, and VarDict).

When provided with a test set, the `classify` subworkflow will produce a plot and metrics to help you evaluate the performance of the trained model and compare it with the other variant callers in the ensemble. These will be stored in the test set's [`prc` folder](#for-test-datasets). The precision-recall plot will be named `results.pdf` and the metrics will be appear in the `pts` directory. You can create the other plots and tables in our paper using the scripts in the [`scripts` directory](/scripts).
