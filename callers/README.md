## Creating a caller script
If you'd like to add another variant caller to the pipeline, you can specify a script to execute it in the [callers dir](https://github.com/aryam7/merge_callers/tree/master/callers).
Caller scripts must have a ".bash" file extension but can be executed using any program by providing the appropriate shebang line. The script should run the caller on a single sample.

Each caller script is referred to in the pipeline by a unique identifier.
Currently, the filename of the script (before the ".bash" extension) is used as the identifier.
When writing your [config file](https://github.com/aryam7/merge_callers/blob/master/config.yaml) you must refer to callers by their identifier.

### Caller script inputs
Each caller script is provided the following arguments respectively.
- a sorted BAM file containing reads for a single sample
- a sorted BED file containing ATAC-seq peaks for a single sample
- a reference genome
- the path to a directory in which the script must place any of its output
- the sample name
- the number of threads to use
- the path to the output directory of a script that prepares files for the caller script (only if it is relevant -- see [below](https://github.com/aryam7/merge_callers/tree/master/callers#caller-scripts-that-depend-on-other-scripts))
- user provided parameters (passed via the [config file](https://github.com/aryam7/merge_callers/blob/master/config.yaml))

### Caller script outputs
Each caller script must use the provided output directory path for all of its output.

Besides this requirement, there is only one other: that at the end of its execution, the script create a TSV file (containing tab-separated values) named by the caller identifier followed by a ".tsv" extension.
The TSV file must include the "CHROM", "POS", "REF", and "ALT" columns first (in that order!).
A header must be present in the TSV file but the "CHROM", "POS", "REF", and "ALT" columns do not need to be named as such.

### Caller scripts that depend on other scripts
Some caller scripts must depend on a different script for special input.
For example, some callers produce files containing both indels and SNVs. The GATK is one such caller.
Separate steps must be performed later to separate the two variant types into different TSV files. How is this achieved?

Well, you can create two different caller scripts `gatk-snp.bash` and `gatk-indel.bash` that each take as input the output of GATK (which is generated using a different special script `gatk.bash`).
Note that `gatk.bash` is not a caller script (and should not be used as one) because it is not expected to produce TSV output.

By providing a dash character `-` in the caller identifier, the caller script (ex: `gatk-snp.bash`) can communicate to the pipeline that it requires input(s) from another special script (ex: `gatk.bash`) with the same filename but without the characters after the final dash (ex: `gatk-snp.bash` => `gatk.bash`).
The pipeline will run this separate script first but with the same parameters as the caller script, and the directory containing its output will be passed to the caller script.

By providing multiple dashes in their caller identifiers using this scheme, the user may design complicated caller script heirarchies involving multiple levels of nesting.
