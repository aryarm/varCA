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
manta_dir="$8"
if [[ -z "$8" ]]; then
	echo "Manta directory not specified. Attempting to retrieve from conda env..." 1>&2
	if [[ -z "$CONDA_PREFIX" ]]; then
		echo "Error: not running in conda env" 1>&2; exit 1;
	else
		manta_dir="$(ls "$CONDA_PREFIX"/share | grep 'manta' | head -n1)"
		if [[ -z "$manta_dir" ]]; then
			echo "Error: coudn't find manta in conda env" 1>&2; exit 1;
		fi
		manta_dir="$CONDA_PREFIX"/share/"$manta_dir"
	fi
else
	manta_dir="$8"
fi
echo "Using manta dir: $manta_dir" 1>&2


# bgzip compress and index the bed file
bgzip -f "$peaks" -c > "$output_dir/peaks.bed.gz" && tabix -p bed -f "$output_dir/peaks.bed.gz"

"$manta_dir"/bin/configManta.py --bam "$bam" --referenceFasta "$genome" --callRegions "$output_dir/peaks.bed.gz" --runDir "$output_dir" --config "$config"
"$output_dir"/runWorkflow.py -m sge -j 36

echo -e "CHROM\tPOS\tREF\tALT\tSVTYPE\tSVLEN\tCIPOS\tCIEND\tHOMLEN\tBND_DEPTH\tMATE_BND_DEPTH\tJUNCTION_QUAL\tGT\tGQ\tPL" > "$output_dir/manta.tsv"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/CIPOS\t%INFO/CIEND\t%INFO/HOMLEN\t%INFO/BND_DEPTH\t%INFO/MATE_BND_DEPTH\t%INFO/JUNCTION_QUAL\t[%GT]\t[%GQ]\t[%PL]\n' "$output_dir"/results/variants/diploidSV.vcf.gz >> "$output_dir/manta.tsv"
