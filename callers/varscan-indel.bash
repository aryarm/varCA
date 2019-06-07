#!/bin/bash

bam="$1"
[[ -z "$1" ]] && { echo "Parameter 1 is empty" ; exit 1; }
peaks="$2"
[[ -z "$2" ]] && { echo "Parameter 2 is empty" ; exit 1; }
genome="$3"
[[ -z "$3" ]] && { echo "Parameter 3 is empty" ; exit 1; }
output="$4"
[[ -z "$4" ]] && { echo "Parameter 4 is empty" ; exit 1; }
output_dir="$5"
[[ -z "$5" ]] && { echo "Parameter 5 is empty" ; exit 1; }

# set the output dir to be non-indel specific
# this allows both snp and indel stuff to be stored in a single "varscan" folder
if [[ "$output_dir" =~ '-indel$' ]]; then
    output_dir="${output_dir%-indel}"
fi

if [ ! -f "$output_dir/varscan.vcf" ]; then
	samtools mpileup -l "$peaks" -f "$genome" -B "$bam" | \
	varscan mpileup2cns --p-value 1 --strand-filter 0 --output-vcf > "$output_dir/varscan.vcf"
fi

echo -e "CHROM\tPOS\tREF\tALT\tADP\tSDP\tDP\tRD\tAD\tFREQ\tPVAL\tABQ\tGT\tGQ" > "$output"
vcftools --vcf "$output_dir/varscan.vcf" --keep-only-indels --recode --recode-INFO-all -c  | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%RC]\t[%CN]\t[%GT]\t[%GQ]\n' - >> "$output"
