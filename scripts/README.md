# scripts
This directory contains various scripts. Some are used by the pipeline.
Others can be used to analyze the tables output by the pipeline.

### configManta.py.ini
The configuration file used by the `manta` variant caller. The path to this file must be specified in the `config.yaml`.

### configureStrelkaGermlineWorkflow.py.ini
The configuration file used by the `strelka` variant caller. The path to this file must also be specified in the `config.yaml`.

### get_cols.bash
A bash script for extracting columns from TSVs via a regex pattern for their column names.

### classify.awk
An awk script for classifying each site in a VCF as DEL, INS, SNP, etc. It accepts a two column table (REF and ALT) from the VCF.
