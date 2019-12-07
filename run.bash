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
#--cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qlog" \
#-j 24 \
#--config output_dir="${out_path}" \
#--latency-wait 60 \
#--use-conda \
#-k \
#"$@" &>"${out_path}/log"

# split the output of the prepare pipeline into a training and testing set
# (only if the output of the prepare pipeline has changed since last time)
#{
#	final="$out_path"/merged_snp/SRR891269
#	([ ! -f "$final"/odd-chrom.tsv.gz ] || [ -n "$(find -L "$final"/final.tsv.gz -prune -newer "$final"/odd-chrom.tsv.gz -exec echo . \;)" ]) && \
#	zcat "$final"/final.tsv.gz | { read -r head && echo "$head" && awk -F'\t' -v 'OFS=\t' '$1 ~ /^[0-9]+$/ && $1%2'; } | gzip > "$final"/odd-chrom.tsv.gz
#	([ ! -f "$final"/even-chrom.tsv.gz ] || [ -n "$(find -L "$final"/final.tsv.gz -prune -newer "$final"/even-chrom.tsv.gz -exec echo . \;)" ]) && \
#	zcat "$out_path"/merged_snp/SRR891269/final.tsv.gz | { read -r head && echo "$head" && awk -F'\t' -v 'OFS=\t' '$1 ~ /^[0-9]+$/ && $1%2 == 0'; } | gzip > "$final"/even-chrom.tsv.gz
#} >>"${out_path}/log" 2>&1
# the above code is only included here so that this file can be executed on the example data all at once
# rather than executing this code, you should manually execute whichever commands will help you prepare the data for the classify pipeline
# (if you're not training the classifier, you probably won't even need to execute any intermediate commands at all)

# classify pipeline -- classify each site; is there a variant there?
snakemake \
-s Snakefiles/Snakefile-classify \
--cluster "qsub -t 1 -V -q iblm.q -j y -o ${out_path}/qlog" \
-j 12 \
--config out="${out_path}/classify" \
--latency-wait 60 \
--use-conda \
-k \
"$@" >>"${out_path}/log" 2>&1

