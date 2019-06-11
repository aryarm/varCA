#!/bin/bash

bam="$1"
[[ -z "$1" ]] && { echo "Parameter 1 is empty" 1>&2; exit 1; }
peaks="$2"
[[ -z "$2" ]] && { echo "Parameter 2 is empty" 1>&2; exit 1; }
genome="$3"
[[ -z "$3" ]] && { echo "Parameter 3 is empty" 1>&2; exit 1; }
output_dir="$4"
[[ -z "$4" ]] && { echo "Parameter 4 is empty" 1>&2; exit 1; }
varscan_dir="$7"
[[ -z "$7" ]] && { echo "Parameter 7 is empty" 1>&2; exit 1; }


echo -e "CHROM\tPOS\tREF\tALT\tADP\tSDP\tDP\tRD\tAD\tFREQ\tPVAL\tABQ\tGT\tGQ" > "$output_dir/varscan-snp.tsv"
vcftools --vcf "$varscan_dir/varscan.vcf" --remove-indels --recode --recode-INFO-all -c | bgzip > "$output_dir/varscan-snp.vcf.gz" && \
tabix -p vcf -f "$output_dir/varscan-snp.vcf.gz"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/ADP\t[%SDP]\t[%DP]\t[%RD]\t[%AD]\t[%FREQ]\t[%PVAL]\t[%ABQ]\t[%GT]\t[%GQ]\n' "$output_dir/varscan-snp.vcf.gz" >> "$output_dir/varscan-snp.tsv"
