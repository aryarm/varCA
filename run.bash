#!/bin/bash
#$ -t 1
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
out_path="$PWD/out"
mkdir -p "$out_path"

# clear leftover log files
if [ -f "${out_path}/log" ]; then
	echo ""> "${out_path}/log";
fi
if [ -f "${out_path}/qlog" ]; then
	echo ""> "${out_path}/qlog";
fi

# make sure that this script is executed from the directory that it lives in!

# Before running these snakemake pipelines, remember to complete both config
# files in the configs/ folder with the required input info. In particular,
# make sure that you have created a samples.tsv file specifying paths to the
# fastq files for each of your samples.

# prepare pipeline -- extract features for each site to prepare for classifying
#snakemake \
#-s Snakefiles/Snakefile-prepare \
#--cluster "qsub -t 1 -V -j y -o ${out_path}/qlog" \
#-j 24 \
#--config output_dir="${out_path}" \
#--latency-wait 60 \
#--use-conda \
#-k \
#"$@" &>"${out_path}/log"

# classify pipeline -- classify each site; is there a variant there?
snakemake \
-s Snakefiles/Snakefile-classify \
--cluster "qsub -t 1 -V -j y -o ${out_path}/qlog" \
-j 12 \
--config out="${out_path}/classify" \
--latency-wait 60 \
--use-conda \
-k \
"$@" >>"${out_path}/log" 2>&1

