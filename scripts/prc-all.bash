#!/bin/bash
#$ -t 1
#$ -q iblm.q
#$ -V
#$ -j y
#$ -cwd

# create all desired precision recall curves for all of our callers

# initialize variables: output dir, input snp table, input indel table, array of snp callers, snp truth caller, array of indel callers, and indel truth caller
out_dir="out"
script_dir="scripts"
snp_table="out/merged_snp/SRR891269.tsv.gz"
IFS=',' read -ra snp_callers <<< "gatk-snp,varscan-snp"
snp_truth="pg-snp"
indel_table="out/merged_indel/SRR891269.tsv.gz"
IFS=',' read -ra indel_callers <<< "gatk-indel,varscan-indel,delly,pindel,manta,strelka"
indel_truth="pg-indel"

for variant_type in 'SNP'; do
	for caller in "${snp_callers[@]}"; do
		"$script_dir"/prc.bash "$snp_table" "$caller" "(REF|ALT)" "$snp_truth" \
		"$variant_type" "$out/prc/points/$caller.${variant_type//,/_}.txt"
	done
done

for variant_type in 'DEL' 'INS' '.' 'DEL,INS'; do
	for caller in "${indel_callers[@]}"; do
		"$script_dir"/prc.bash "$indel_table" "$caller" "(REF|ALT)" "$indel_truth" \
		"$variant_type" "$out/prc/points/$caller.${variant_type//,/_}.txt"
	done
done
