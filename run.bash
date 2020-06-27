#!/bin/bash
#$ -t 1
#$ -V
#$ -j y
#$ -cwd
#$ -o /dev/null
#$ -e /dev/null

# An example bash script demonstrating how to run the entire snakemake pipeline
# This script creates two separate log files:
# 	1) log - the basic snakemake log of completed rules
# 	2) qlog - a more detailed log of the progress of each rule and any errors

# Before running the snakemake pipeline, remember to complete the config.yaml
# file in the configs/ folder with the required input info. In particular,
# make sure that you have created a samples.tsv file specifying paths to the
# fastq (or bam) files for each of your samples.
# Make sure that this script is executed from the directory that it lives in!

out_path="out" # you can specify a dir for all output here (or in the configs)
mkdir -p "$out_path"

# clear leftover log files
if [ -f "$out_path/log" ]; then
	echo ""> "$out_path/log";
fi
if [ -f "$out_path/qlog" ]; then
	echo ""> "$out_path/qlog";
fi

# check: are we being executed from within qsub?
if [ "$ENVIRONMENT" = "BATCH" ]; then
	snakemake \
	--cluster "qsub -t 1 -V -j y -cwd -o $out_path/qlog" \
	--config out="$out_path" \
	--latency-wait 60 \
	--use-conda \
	-k \
	-j \
	"$@" &>"$out_path/log"

	# -----------
	# An example running the classify subworkflow
	# Remember to complete the classify.yaml config file before running this!
	# snakemake \
	# -s rules/classify.smk \
	# --cluster "qsub -t 1 -V -j y -cwd -o $out_path/qlog" \
	# --config out="$out_path/classify" \
	# --latency-wait 60 \
	# --use-conda \
	# -k \
	# -j \
	# "$@" &>"$out_path/log"
else
	snakemake \
	--config out="$out_path" \
	--latency-wait 60 \
	--use-conda \
	-k \
	-j \
	"$@" 2>"$out_path/log" >"$out_path/qlog"

	# -----------
	# An example running the classify subworkflow
	# Remember to complete the classify.yaml config file before running this!
	# snakemake \
	# -s rules/classify.smk \
	# --config out="$out_path/classify" \
	# --latency-wait 60 \
	# --use-conda \
	# -k \
	# -j \
	# "$@" 2>"$out_path/log" >"$out_path/qlog"
fi

