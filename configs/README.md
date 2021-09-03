## Completing the configuration files
Before running the pipeline or its subworkflows, you must complete the relevant confiugration files to specify inputs, set parameters, and indicate where to dump output. Read about each config file below to learn about whether and how you should fill it out.
You can override any of the entries in the config files at execution time by passing them via the command line. For example, to set the output directory of either pipeline, specify `--config out="${out_path}"` as a parameter to Snakemake.

### [config.yaml](config.yaml)
Edit the configuration in this file if you want to execute the master pipeline on your own data. It is currently configured to run on the example data provided.

You must provide the main inputs in this config file:

1. Paired ATAC-seq reads as FASTQ (or BAM) files, specified via a tab-delimited "samples file"
2. A properly indexed reference genome for the samples in #1
3. Lists of the variant callers to use in the ensemble when predicting a) SNVs and b) indels
    - The recommended callers have already been filled out for you
    - Each caller must be referred to by its caller ID (see the [caller README](/callers/README.md) for more information)
    - You should provide these callers in a specific order, such that the callers that are more likely to make an accurate call are listed first
4. Trained, classification models that can be used for predicting a) SNVs and b) indels from the VCF output of each variant caller
    - Currently, the trained models provided with the example data are used. These have been trained on Platinum Genomes from GM12878 (a high quality truth set), but you can also [create your own trained models](/rules#creating-your-own-trained-model).

If you provide BAM files instead of FASTQs, the pipeline will assume that you have already removed PCR duplicates, so it will skip the duplicate removal step. If you would like to skip the peak calling step as well, you will also need to provide a BED file (with a .bed extension) containing the peaks that should be used. Note also that there are many requirements that your BAM files must adhere to:

1. They must have [read group information](https://gatk.broadinstitute.org/hc/articles/360035890671)
2. They cannot have PCR duplicates
3. They must have '.bam' file extensions
4. There must be an index ('.bam.bai' file) in the same directory

The reference genome must be properly indexed. The pipeline needs [BWA index files](https://hcc.unl.edu/docs/applications/app_specific/bioinformatics_tools/alignment_tools/bwa/running_bwa_commands/#bwa-index), a [GATK `.dict` file](https://gatk.broadinstitute.org/hc/en-us/articles/360037068312-CreateSequenceDictionary-Picard-), and [a samtools index file](http://www.htslib.org/doc/samtools-faidx.html).

You may also provide any other configuration options from the `prepare.yaml` config file.

Once you have finished filling out the `config.yaml` file, you should be ready to execute the master pipeline on your own data. Just execute [the `run.bash` script](/run.bash).

The output of the master pipeline is the `final.vcf.gz` file in the `classify` directory. For example, when the master pipeline is executed on the example data, the output is `out/classify/*/final.vcf.gz`.

### [prepare.yaml](prepare.yaml)
The `prepare` subworkflow executes each variant caller in the ensemble and prepares the combined output for use by the `classify` subworkflow.

The master pipeline automatically executes the `prepare` subworkflow using the configuration in `config.yaml`. However, it is also possible to execute the `prepare` subworkflow on its own, independent of the master pipeline. See the [rules README](/rules/README.md) for instructions on how to do this. You must first complete the configuration in `prepare.yaml`. At the moment, it is partially configured for the GM12878 data from [Buenrestro et al](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753).

The main inputs in this config file are:

1. Paired ATAC-seq reads as FASTQ (or BAM) files, specified via a tab-delimited "samples file"
2. A properly indexed reference genome for the samples in #1
3. Lists of the variant callers to use in the ensemble when predicting a) SNVs and b) indels
    - Each caller must be referred to by its caller ID (see the [caller README](/callers/README.md) for more information)

As described for `config.yaml` above, you can provide a BED file if you want to skip the peak calling step. All of the requirements listed for the BAM files in `config.yaml` also apply here.

You may also provide the name of a directory in which to store all output. Any other configuration options not listed here have reasonable defaults and are described in the [`prepare.yaml` config file](prepare.yaml).

The primary output of the subworkflow is a gzipped TSV containing VCF columns from every specified variant caller.

You can choose which variant callers are used by providing their "caller IDs" in `prepare.yaml`. If you'd like to include a caller that we haven't implemented yet, simply follow these steps:

1. Create [a caller script](/callers#creating-a-caller-script) that describes how the caller should be executed
2. Add its caller ID to the config files
3. Provide any [special parameters](callers.yaml) to your caller script or specify any special VCF column names that you'd like to train the classifier on
4. Create new, trained classification models using the training data (see the [rules README](/rules#creating-your-own-trained-model) for instructions)

### [callers.yaml](callers.yaml)
The `prepare` subworkflow (and, by extension, the master pipeline) read parameters specific to each variant caller from the `callers.yaml` config file. These parameters change the data that is output by the `prepare` subworkflow. As a general rule, changing any of the config options in the `callers.yaml` file will require that you [retrain the classification models](/rules#training-and-testing-varca) to recognize the new data that you are providing to the `classify` subworkflow.

For [each caller script in the callers directory](/callers) you can specify:

1. The names of the fields in the VCFs generated by each caller script which should be included in the final TSV output by the `prepare` subworkflow
    - The final TSV is used for training, testing, and applying the classifier. We recommend using fields that are numerical rather than categorical (although both are acceptable).
2. Any extra parameters that should be passed to the caller scripts when they are executed
    - For example, some callers (Strelka, Manta, and BreakCA) require that you provide the path to executables on your system if you aren't running snakemake with the `--use-conda` parameter
3. The type of output that the caller script generates (if not a VCF)
4. If a caller script outputs N/A values for some sites, what should those N/A values be replaced with, if not 0?

### [classify.yaml](classify.yaml)
The `classify` subworkflow can be used to train, test, and apply models capable of classifying sites by their variant type (ie SNP, DEL, INS, no variant, etc). It uses datasets prepared by the `prepare` subworkflow.

The master pipeline automatically executes the `classify` subworkflow using the configuration in `config.yaml`. However, it is also possible to execute the `classify` subworkflow on its own, independent of the master pipeline. See the [rules README](/rules/README.md) for instructions on how to do this. You must first complete the configuration in `classify.yaml`. At the moment, it is partially configured for the GM12878 data from [Buenrestro et al](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47753).

The main inputs in this config file are:

1. The TSV datasets created by the `prepare` pipeline
    - Unless a trained classification model is already provided, one of these will be used to create the model, while the rest will be used for either testing or applying the trained model. You can specify which should be used for which in the "train" and "predict" config options.
2. A properly indexed reference genome for the samples in #1
3. A pre-trained classification model, if the training step should be skipped

You can also provide the name of a directory in which to store all output. Any other configuration options not listed here have reasonable defaults and are described in the [`classify.yaml` config file](classify.yaml).

The primary output of the subworkflow is a gzipped VCF containing recalibrated QUAL scores based on the probabilities output by the classifier. The VCF is constructed by sampling variants from the variant callers in the ensemble, in an order that you can specify in `classify.yaml`. For each variant, we indicate which caller was used as a field in INFO.

## Other config files

### [configManta.py.ini](configManta.py.ini)
This is the configuration file used by the `manta` variant caller. The path to this file must be specified in `callers.yaml`.

### [configureStrelkaGermlineWorkflow.py.ini](configureStrelkaGermlineWorkflow.py.ini)
This is the configuration file used by the `strelka` variant caller. The path to this file must also be specified in `callers.yaml`.
