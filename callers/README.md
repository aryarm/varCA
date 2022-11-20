## What is a caller script?
An isolated script in which you can define the logic for running a variant caller. These scripts are called, in turn, by the `prepare` subworkflow.

## Creating a caller script
If you'd like to add another variant caller to the pipeline, you can specify a script to execute it in the [callers dir](/callers).
The script should run the caller on a single sample.

Each caller script is referred to in the pipeline by a unique identifier.
Currently, the filename of the script is used as the identifier.
When writing the `prepare.yaml` and `callers.yaml` [config files](/configs) for the `prepare` subworkflow, you must refer to callers by their identifier.

### Caller script inputs
Each caller script is provided the following arguments respectively (any of which can be ignored):
- a sorted BAM file containing reads for a single sample
- a sorted BED file containing ATAC-seq peaks for a single sample
- a reference genome
- the path to a directory in which the script must place any of its output
- the sample name
- the path to the output directory of a script that prepares files for the caller script (only if it is relevant -- see [below](#caller-scripts-that-depend-on-other-scripts))
- user provided parameters (passed via the [`callers.yaml` config file](/configs/callers.yaml))

### Caller script outputs
Each caller script must use the provided output directory path for all of its output.

Besides this requirement, there is only one other: that at the end of its execution, the script create a VCF file named by the caller identifier followed by a ".vcf" extension. The VCF file cannot be gzipped.

You may optionally specify that the caller script outputs a TSV instead of a VCF. Just add an 'ext' attribute with a value of 'tsv' to the caller-specific parameters in the [`callers.yaml` config file](/configs/callers.yaml). The TSV file must be named similarly to the VCF, except that it must have a ".tsv" instead of a ".vcf" extension. The first two columns of the TSV must be "CHROM" and "POS", in that order.

### Caller scripts that depend on other scripts
Some caller scripts must depend on a different script for special input.

For example, some callers produce files containing both indels and SNVs. The GATK is one such caller.
Separate steps must be performed later to separate the two variant types into different TSV files. How is this achieved?

Well, you can create two different caller scripts `gatk-snp` and `gatk-indel` that each take as input the output of GATK (which is generated using a different special script `gatk`).
Then `gatk-snp` and `gatk-indel` can each extract SNVs and indels from the output of `gatk`, respectively.
Note that in this example, `gatk` is not a caller script (and cannot be used as one) because it will not produce the correct output. Even if it did, it would be ignored.

By providing a dash character `-` in the caller identifier, the caller script (ex: `gatk-snp`) can communicate to the pipeline that it requires input(s) from another special script with the same filename but without the characters after the final dash (ex: `gatk-snp` => `gatk`).
The pipeline will run this separate script first but with the same parameters as the caller script, and the directory containing its output will be passed as an argument to the original caller script. (For example, the directory containing the output of the `gatk` special script will be passed as a parameter to the `gatk-snp` caller script.)

By providing multiple dashes in your caller identifiers using this scheme, you may design complicated caller script heirarchies involving multiple levels of nesting. In a sense, you can create small pipelines within the pipeline.

### Caller script dependencies
If you caller script requires specific dependencies, just add them to the [prepare.yml](/envs/prepare.yml) [conda environment file](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually).
Any packages listed in this file will be available to your caller script when it is executed.

If you would like to use a newer version of an existing variant caller, just change the version specified in the [prepare.yml](/envs/prepare.yml) file. Snakemake will automatically update the package upon your next execution of the pipeline.
