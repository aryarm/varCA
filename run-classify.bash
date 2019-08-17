#!/bin/bash
#$ -t 1
#$ -q iblm.q
#$ -V
#$ -j y
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null

# An example bash script demonstrating how to run the entire snakemake pipeline
# on an SGE cluster
# This script creates two separate log files:
# 	1) log - the basic snakemake log of completed rules
# 	2) qlog - a more detailed log of the progress of each rule and any errors

# you can specify a directory for all output here:
out_path="$PWD/out/classify"
mkdir -p "$out_path"

# clear leftover log files
if [ -f "${out_path}/log" ]; then
	echo ""> "${out_path}/log";
fi
if [ -f "${out_path}/qlog" ]; then
	echo ""> "${out_path}/qlog";
fi

# make sure that this script is executed from the directory that it lives in!

# Before running this snakemake pipeline, remember to verify that the config.yaml
# file has been appropriately completed with the required input info. In
# particular, make sure that you have created a samples.tsv file specifying
# paths to the fastq files for each of your samples.

snakemake \
-s Snakefiles/Snakefile-classify \
--cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qlog" \
-j 12 \
--config out="${out_path}" \
--latency-wait 60 \
--use-conda \
-k \
--notemp \
"$@" &>"${out_path}/log"
