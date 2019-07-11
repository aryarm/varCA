#!/bin/bash
#$ -t 1
#$ -q iblm.q
#$ -V
#$ -j y
#$ -cwd

# create all desired precision recall curves for all of our callers

# initialize variables: output dir, input snp table, input indel table, array of snp callers, snp truth caller, array of indel callers, and indel truth caller
out="out"
script_dir="scripts"
snp_table="out/merged_snp/SRR891269.tsv.gz"
IFS=',' read -ra snp_callers <<< "gatk-snp,varscan-snp,vardict-snp"
snp_truth="pg-snp"
indel_table="out/merged_indel/SRR891269.tsv.gz"
IFS=',' read -ra indel_callers <<< "gatk-indel,varscan-indel,vardict-indel,pindel,illumina-manta,illumina-strelka"
indel_truth="pg-indel"

mkdir -p "$out/prc/points"
mkdir -p "$out/prc/points/curves"
mkdir -p "$out/prc/points/singles"
mkdir -p "$out/prc/plots"

for variant_type in '.'; do
	for caller in "${snp_callers[@]}"; do
		flip='0'
		sort='1'
		if [ "$caller" == 'gatk-snp' ]; then
			col="QD"
		elif [ "$caller" == "varscan-snp" ]; then
			col="PVAL"
			flip='1'
		elif [ "$caller" == "vardict-snp" ]; then
			col="QUAL"
		else
			col="QG"
		fi
		[ ! -f "$out/prc/points/curves/$caller.${variant_type//,/_}.txt" ] && \
		echo "$caller: prc.bash" 1>&2 && \
		"$script_dir"/prc.bash "$snp_table" "$caller~$col" "$snp_truth" \
		"$variant_type" "$out/prc/points/curves/$caller.${variant_type//,/_}.txt" 'gatk-snp~DP>10' "$flip" "$sort"
		[ ! -f "$out/prc/points/singles/$caller.${variant_type//,/_}.txt" ] && \
		echo "$caller: metrics.bash" 1>&2 && \
		"$script_dir"/metrics.bash "$snp_table" "$caller" "$snp_truth" \
		"$variant_type" "$out/prc/points/singles/$caller.${variant_type//,/_}.txt" 'gatk-snp~DP>10'
	done
done

# for variant_type in '.' 'DEL' 'INS' 'DEL,INS'; do
for variant_type in '.'; do
	for caller in "${indel_callers[@]}"; do
		flip='0'
		sort='1'
		if [ "$caller" == 'gatk-indel' ]; then
			col="QD"
		elif [ "$caller" == "varscan-indel" ]; then
			col="PVAL"
			flip='1'
		elif [ "$caller" == "vardict-indel" ]; then
			col="QUAL"
		elif [ "$caller" == 'pindel' ] || [ "$caller" == 'illumina-manta' ] || [ "$caller" == 'illumina-strelka' ]; then
			col=''
		else
			col="QG"
		fi
		if [ ! -z "$col" ]; then
			[ ! -f "$out/prc/points/curves/$caller.${variant_type//,/_}.txt" ] && \
			echo "$caller: prc.bash" 1>&2 && \
			"$script_dir"/prc.bash "$indel_table" "$caller~$col" "$indel_truth" \
			"$variant_type" "$out/prc/points/curves/$caller.${variant_type//,/_}.txt" 'gatk-indel~DP>10' "$flip" "$sort"
		fi
		[ ! -f "$out/prc/points/singles/$caller.${variant_type//,/_}.txt" ] && \
		echo "$caller: metrics.bash" 1>&2 && \
		"$script_dir"/metrics.bash "$indel_table" "$caller" "$indel_truth" \
		"$variant_type" "$out/prc/points/singles/$caller.${variant_type//,/_}.txt" 'gatk-indel~DP>10'
	done
done

[ ! -f "$out/prc/plots/all-callers.snv.png" ] && \
python "$script_dir"/prc.py "$out/prc/plots/all-callers.snv.png" \
--gatk-snp "$out/prc/points/curves/gatk-snp...txt" \
--gatk-snp-pt "$out/prc/points/singles/gatk-snp...txt" \
--varscan-snp "$out/prc/points/curves/varscan-snp...txt" \
--varscan-snp-pt "$out/prc/points/singles/varscan-snp...txt" \
--vardict-snp "$out/prc/points/curves/vardict-snp...txt" \
--vardict-snp-pt "$out/prc/points/singles/vardict-snp...txt"

[ ! -f "$out/prc/plots/all-callers.indel.png" ] && \
python "$script_dir"/prc.py "$out/prc/plots/all-callers.indel.png" \
--gatk-indel "$out/prc/points/curves/gatk-indel...txt" \
--gatk-indel-pt "$out/prc/points/singles/gatk-indel...txt" \
--varscan-indel "$out/prc/points/curves/varscan-indel...txt" \
--varscan-indel-pt "$out/prc/points/singles/varscan-indel...txt" \
--vardict-indel "$out/prc/points/curves/vardict-indel...txt" \
--vardict-indel-pt "$out/prc/points/singles/vardict-indel...txt" \
--pindel-pt "$out/prc/points/singles/pindel...txt" \
--illumina-manta-pt "$out/prc/points/singles/illumina-manta...txt" \
--illumina-strelka-pt "$out/prc/points/singles/illumina-strelka...txt"

# [ ! -f "$out/prc/plots/all-callers-micro.indel.png" ] && \
# python "$script_dir"/prc.py "$out/prc/plots/all-callers-micro.indel.png" \
# --gatk-indel "$out/prc/points/curves/gatk-indel.DEL_INS.txt" \
# --gatk-indel-pt "$out/prc/points/singles/gatk-indel.DEL_INS.txt" \
# --varscan-indel "$out/prc/points/curves/varscan-indel.DEL_INS.txt" \
# --varscan-indel-pt "$out/prc/points/singles/varscan-indel.DEL_INS.txt" \
# --vardict-indel "$out/prc/points/curves/vardict-indel.DEL_INS.txt" \
# --vardict-indel-pt "$out/prc/points/singles/vardict-indel.DEL_INS.txt" \
# --delly "$out/prc/points/curves/delly.DEL_INS.txt" \
# --delly-pt "$out/prc/points/singles/delly.DEL_INS.txt" \
# --pindel-pt "$out/prc/points/singles/pindel.DEL_INS.txt" \
# --illumina-manta "$out/prc/points/curves/illumina-manta.DEL_INS.txt" \
# --illumina-manta-pt "$out/prc/points/singles/illumina-manta.DEL_INS.txt" \
# --illumina-strelka "$out/prc/points/curves/illumina-strelka.DEL_INS.txt" \
# --illumina-strelka-pt "$out/prc/points/singles/illumina-strelka.DEL_INS.txt"

# [ ! -f "$out/prc/plots/all-callers-micro.indel.png" ] && \
# python "$script_dir"/prc.py "$out/prc/plots/all-callers-micro.indel.png" \
# --gatk-indel "$out/prc/points/curves/gatk-indel.DEL_INS.txt" \
# --varscan-indel "$out/prc/points/curves/varscan-indel.DEL_INS.txt" \
# --vardict-indel "$out/prc/points/curves/vardict-indel.DEL_INS.txt" \
# --illumina-manta "$out/prc/points/curves/illumina-manta.DEL_INS.txt" \
# --illumina-strelka "$out/prc/points/curves/illumina-strelka.DEL_INS.txt"

# [ ! -f "$out/prc/plots/gatk-indel-multi.png" ] && \
# python "$script_dir"/prc.py "$out/prc/plots/gatk-indel-multi.png" \
# --gatk-micro "$out/prc/points/curves/gatk-indel.DEL_INS.txt" \
# --gatk-ins "$out/prc/points/curves/gatk-indel.INS.txt" \
# --gatk-del "$out/prc/points/curves/gatk-indel.DEL.txt"

# [ ! -f "$out/prc/plots/varscan-indel-multi.png" ] && \
# python "$script_dir"/prc.py "$out/prc/plots/varscan-indel-multi.png" \
# --varscan-micro "$out/prc/points/curves/varscan-indel.DEL_INS.txt" \
# --varscan-ins "$out/prc/points/curves/varscan-indel.INS.txt" \
# --varscan-del "$out/prc/points/curves/varscan-indel.DEL.txt"

# [ ! -f "$out/prc/plots/vardict-indel-multi.png" ] && \
# python "$script_dir"/prc.py "$out/prc/plots/vardict-indel-multi.png" \
# --vardict-micro "$out/prc/points/curves/vardict-indel.DEL_INS.txt" \
# --vardict-ins "$out/prc/points/curves/vardict-indel.INS.txt" \
# --vardict-del "$out/prc/points/curves/vardict-indel.DEL.txt"
