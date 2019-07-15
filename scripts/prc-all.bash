#!/bin/bash
#$ -t 1
#$ -q iblm.q
#$ -V
#$ -j y
#$ -cwd

# create all desired precision recall curves for all of our callers

# initialize variables: output dir, input snp table, input indel table, array of snp callers, snp truth caller, array of indel callers, and indel truth caller
out="out/prc"
script_dir="scripts"
snp_table="out/merged_snp/SRR891269.tsv.gz"
IFS=',' read -ra snp_callers <<< "gatk-snp,varscan-snp,vardict-snp"
snp_truth="pg-snp"
indel_table="out/merged_indel/SRR891269.tsv.gz"
IFS=',' read -ra indel_callers <<< "gatk-indel,varscan-indel,vardict-indel,pindel,illumina-manta,illumina-strelka"
indel_truth="pg-indel"

# what depths would we like to filter at?
depths=('' 0 5 10 20)

mkdir -p "$out/points"
mkdir -p "$out/plots"

for depth in "${depths[@]}"; do
	if [ -z "$depth" ]; then
		depth_expr=''
	else
		depth_expr='gatk-snp~DP>'"$depth"
	fi
	out_dir="$out/points/d$depth"
	mkdir -p "$out_dir"
	mkdir -p "$out_dir/curves"
	mkdir -p "$out_dir/singles"
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
			[ ! -f "$out_dir/curves/$caller.${variant_type//,/_}.txt" ] && \
			echo "$caller: prc.bash at d=$depth" 1>&2 && \
			"$script_dir"/prc.bash "$snp_table" "$caller~$col" "$snp_truth" \
			"$variant_type" "$out_dir/curves/$caller.${variant_type//,/_}.txt" "$depth_expr" "$flip" "$sort"
			[ ! -f "$out_dir/singles/$caller.${variant_type//,/_}.txt" ] && \
			echo "$caller: metrics.bash at d=$depth" 1>&2 && \
			"$script_dir"/metrics.bash "$snp_table" "$caller" "$snp_truth" \
			"$variant_type" "$out_dir/singles/$caller.${variant_type//,/_}.txt" "$depth_expr"
		done
		# create a summary table of all single point metrics
		if [ "$variant_type" == "." ]; then
			single_paths=( "${snp_callers[@]/#/$out_dir/singles/}" )
			[ ! -f "$out_dir/singles/snp_summary.tsv" ] && \
			echo "snp_summary at d=$depth" 1>&2 && \
			echo -e "recall\nprecision\nf-beta\npositives\nnegatives" | paste - "${single_paths[@]/%/...txt}" > "$out_dir/singles/snp_summary.tsv" && \
			echo "$(echo 'metric' "${snp_callers[@]}" | tr ' ' '\t' | cat - "$out_dir/singles/snp_summary.tsv")" > "$out_dir/singles/snp_summary.tsv"
		fi
	done
done

for depth in "${depths[@]}"; do
	if [ -z "$depth" ]; then
		depth_expr=''
		variant_types=('.' 'DEL' 'INS' 'DEL,INS')
	else
		depth_expr='gatk-indel~DP>'"$depth"
		variant_types=('.')
	fi
	out_dir="$out/points/d$depth"
	mkdir -p "$out_dir"
	mkdir -p "$out_dir/curves"
	mkdir -p "$out_dir/singles"
	for variant_type in "${variant_types[@]}"; do
		for caller in "${indel_callers[@]}"; do
			flip='0'
			sort='1'
			if [ "$caller" == 'gatk-indel' ]; then
				col="QD"
			elif [ "$caller" == "varscan-indel" ]; then
				col="PVAL"
				flip='1'
			elif [ "$caller" == "vardict-indel" ] || [ "$caller" == 'illumina-manta' ] || [ "$caller" == 'illumina-strelka' ]; then
				col="QUAL"
			elif [ "$caller" == 'pindel' ]; then
				col=''
			else
				col="QG"
			fi
			if [ ! -z "$col" ]; then
				[ ! -f "$out_dir/curves/$caller.${variant_type//,/_}.txt" ] && \
				echo "$caller: prc.bash at d=$depth" 1>&2 && \
				"$script_dir"/prc.bash "$indel_table" "$caller~$col" "$indel_truth" \
				"$variant_type" "$out_dir/curves/$caller.${variant_type//,/_}.txt" "$depth_expr" "$flip" "$sort"
			fi
			[ ! -f "$out_dir/singles/$caller.${variant_type//,/_}.txt" ] && \
			echo "$caller: metrics.bash at d=$depth" 1>&2 && \
			"$script_dir"/metrics.bash "$indel_table" "$caller" "$indel_truth" \
			"$variant_type" "$out_dir/singles/$caller.${variant_type//,/_}.txt" "$depth_expr"
		done
		# create a summary table of all single point metrics
		if [ "$variant_type" == "." ]; then
			single_paths=( "${indel_callers[@]/#/$out_dir/singles/}" )
			[ ! -f "$out_dir/singles/indel_summary.tsv" ] && \
			echo "indel_summary at d=$depth" 1>&2 && \
			echo -e "recall\nprecision\nf-beta\npositives\nnegatives" | paste - "${single_paths[@]/%/...txt}" > "$out_dir/singles/indel_summary.tsv" && \
			echo "$(echo 'metric' "${indel_callers[@]}" | tr ' ' '\t' | cat - "$out_dir/singles/indel_summary.tsv")" > "$out_dir/singles/indel_summary.tsv"
			if [ -f "$out_dir/singles/snp_summary.tsv" ] && [ -f "$out_dir/singles/indel_summary.tsv" ] && [ ! -f "$out_dir/singles/summary.tsv" ]; then
				echo "overall summary at d=$depth" 1>&2
				cut -f 2- "$out_dir/singles/indel_summary.tsv" | paste "$out_dir/singles/snp_summary.tsv" - > "$out_dir/singles/summary.tsv"
			fi
		fi
	done
done


for depth in "${depths[@]}"; do
	points_out="$out/points/d$depth"
	plots_out="$out/plots/d$depth"
	mkdir -p "$plots_out"

	[ ! -f "$plots_out/all-callers.snv.png" ] && \
	echo "all-callers.snv at d=$depth" 1>&2 && \
	python "$script_dir"/prc.py "$plots_out/all-callers.snv.png" \
	--gatk-snp "$points_out/curves/gatk-snp...txt" \
	--gatk-snp-pt "$points_out/singles/gatk-snp...txt" \
	--varscan-snp "$points_out/curves/varscan-snp...txt" \
	--varscan-snp-pt "$points_out/singles/varscan-snp...txt" \
	--vardict-snp "$points_out/curves/vardict-snp...txt" \
	--vardict-snp-pt "$points_out/singles/vardict-snp...txt"

	[ ! -f "$plots_out/all-callers.indel.png" ] && \
	echo "all-callers.indel at d=$depth" 1>&2 && \
	python "$script_dir"/prc.py "$plots_out/all-callers.indel.png" \
	--gatk-indel "$points_out/curves/gatk-indel...txt" \
	--gatk-indel-pt "$points_out/singles/gatk-indel...txt" \
	--varscan-indel "$points_out/curves/varscan-indel...txt" \
	--varscan-indel-pt "$points_out/singles/varscan-indel...txt" \
	--vardict-indel "$points_out/curves/vardict-indel...txt" \
	--vardict-indel-pt "$points_out/singles/vardict-indel...txt" \
	--pindel-pt "$points_out/singles/pindel...txt" \
	--illumina-manta "$points_out/curves/illumina-manta...txt" \
	--illumina-manta-pt "$points_out/singles/illumina-manta...txt" \
	--illumina-strelka "$points_out/curves/illumina-strelka...txt" \
	--illumina-strelka-pt "$points_out/singles/illumina-strelka...txt"

	if [ -z "$depth" ]; then
		[ ! -f "$plots_out/all-callers-micro.indel.png" ] && \
		echo "all-callers-micro.indel" 1>&2 && \
		python "$script_dir"/prc.py "$plots_out/all-callers-micro.indel.png" \
		--gatk-indel "$points_out/curves/gatk-indel.DEL_INS.txt" \
		--gatk-indel-pt "$points_out/singles/gatk-indel.DEL_INS.txt" \
		--varscan-indel "$points_out/curves/varscan-indel.DEL_INS.txt" \
		--varscan-indel-pt "$points_out/singles/varscan-indel.DEL_INS.txt" \
		--vardict-indel "$points_out/curves/vardict-indel.DEL_INS.txt" \
		--vardict-indel-pt "$points_out/singles/vardict-indel.DEL_INS.txt" \
		--pindel-pt "$points_out/singles/pindel.DEL_INS.txt" \
		--illumina-manta "$points_out/curves/illumina-manta.DEL_INS.txt" \
		--illumina-manta-pt "$points_out/singles/illumina-manta.DEL_INS.txt" \
		--illumina-strelka "$points_out/curves/illumina-strelka.DEL_INS.txt" \
		--illumina-strelka-pt "$points_out/singles/illumina-strelka.DEL_INS.txt"

		[ ! -f "$plots_out/all-callers-micro.indel.png" ] && \
		echo "all-callers-micro.indel" 1>&2 && \
		python "$script_dir"/prc.py "$plots_out/all-callers-micro.indel.png" \
		--gatk-indel "$points_out/curves/gatk-indel.DEL_INS.txt" \
		--varscan-indel "$points_out/curves/varscan-indel.DEL_INS.txt" \
		--vardict-indel "$points_out/curves/vardict-indel.DEL_INS.txt" \
		--illumina-manta "$points_out/curves/illumina-manta.DEL_INS.txt" \
		--illumina-strelka "$points_out/curves/illumina-strelka.DEL_INS.txt"

		[ ! -f "$plots_out/gatk-indel-multi.png" ] && \
		echo "gatk-indel-multi" 1>&2 && \
		python "$script_dir"/prc.py "$plots_out/gatk-indel-multi.png" \
		--gatk-micro "$points_out/curves/gatk-indel.DEL_INS.txt" \
		--gatk-ins "$points_out/curves/gatk-indel.INS.txt" \
		--gatk-del "$points_out/curves/gatk-indel.DEL.txt"

		[ ! -f "$plots_out/varscan-indel-multi.png" ] && \
		echo "varscan-indel-multi" 1>&2 && \
		python "$script_dir"/prc.py "$plots_out/varscan-indel-multi.png" \
		--varscan-micro "$points_out/curves/varscan-indel.DEL_INS.txt" \
		--varscan-ins "$points_out/curves/varscan-indel.INS.txt" \
		--varscan-del "$points_out/curves/varscan-indel.DEL.txt"

		[ ! -f "$plots_out/vardict-indel-multi.png" ] && \
		echo "vardict-indel-multi" 1>&2 && \
		python "$script_dir"/prc.py "$plots_out/vardict-indel-multi.png" \
		--vardict-micro "$points_out/curves/vardict-indel.DEL_INS.txt" \
		--vardict-ins "$points_out/curves/vardict-indel.INS.txt" \
		--vardict-del "$points_out/curves/vardict-indel.DEL.txt"

		[ ! -f "$plots_out/illumina-manta-multi.png" ] && \
		echo "illumina-manta-multi" 1>&2 && \
		python "$script_dir"/prc.py "$plots_out/illumina-manta-multi.png" \
		--manta-micro "$points_out/curves/illumina-manta.DEL_INS.txt" \
		--manta-ins "$points_out/curves/illumina-manta.INS.txt" \
		--manta-del "$points_out/curves/illumina-manta.DEL.txt"

		[ ! -f "$plots_out/illumina-strelka-multi.png" ] && \
		echo "illumina-strelka-multi" 1>&2 && \
		python "$script_dir"/prc.py "$plots_out/illumina-strelka-multi.png" \
		--strelka-micro "$points_out/curves/illumina-strelka.DEL_INS.txt" \
		--strelka-ins "$points_out/curves/illumina-strelka.INS.txt" \
		--strelka-del "$points_out/curves/illumina-strelka.DEL.txt"
	fi
done
