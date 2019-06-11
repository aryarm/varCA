#!/bin/bash


bam="$1"
[[ -z "$1" ]] && { echo "Parameter 1 is empty" 1>&2; exit 1; }
peaks="$2"
[[ -z "$2" ]] && { echo "Parameter 2 is empty" 1>&2; exit 1; }
genome="$3"
[[ -z "$3" ]] && { echo "Parameter 3 is empty" 1>&2; exit 1; }
output_dir="$4"
[[ -z "$4" ]] && { echo "Parameter 4 is empty" 1>&2; exit 1; }
config="$7"
[[ -z "$7" ]] && { echo "Parameter 7 is empty" 1>&2; exit 1; }
strelka_dir="$8"
if [[ -z "$8" ]]; then
	echo "Strelka directory not specified. Attempting to retrieve from conda env..." 1>&2
	if [[ -z "$CONDA_PREFIX" ]]; then
		echo "Error: not running in conda env" 1>&2; exit 1;
	else
		strelka_dir="$(ls "$CONDA_PREFIX"/share | grep 'strelka' | head -n1)"
		if [[ -z "$strelka_dir" ]]; then
			echo "Error: coudn't find strelka in conda env" 1>&2; exit 1;
		fi
		strelka_dir="$CONDA_PREFIX"/share/"$strelka_dir"
	fi
else
	strelka_dir="$8"
fi
echo "Using strelka dir: $strelka_dir" 1>&2


# bgzip compress and index the bed file
bgzip -f "$peaks" -c > "$output_dir/peaks.bed.gz" && tabix -p bed -f "$output_dir/peaks.bed.gz"

"$strelka_dir"/bin/configureStrelkaGermlineWorkflow.py --bam "$bam" --referenceFasta "$genome" --callRegions "$output_dir/peaks.bed.gz" --runDir "$output_dir" --config "$config"
"$output_dir"/runWorkflow.py -m sge -j 36

echo -e "CHROM\tPOS\tREF\tALT\tSNVHPOL\tRU\tREFREP\tIDREP\tMQ\tGT\tGQ\tGQX\tDP\tDPF\tMIN_DP\tAD\tADF\tADR\tDPI\tPL\tPS\tSB" > "$output_dir/strelka.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SNVHPOL\t%INFO/RU\t%INFO/REFREP\t%INFO/IDREP\t%INFO/MQ\t[%GT]\t[%GQ\t[%GQX]\t[%DP]\t[%DPF]\t[%MIN_DP]\t[%AD]\t[%ADF]\t[%ADR]\t[%DPI]\t[%PL]\t[%PS]\t[%SB]\n' "$output_dir"/results/variants/genome.S1.vcf.gz >> "$output_dir/strelka.tsv"
