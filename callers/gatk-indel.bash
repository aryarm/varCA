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
# this allows both snp and indel stuff to be stored in a single "gatk" folder
if [[ "$output_dir" =~ '-indel$' ]]; then
    output_dir="${output_dir%-indel}"
fi

# run gatk haplotype caller if needed
if [ ! -f "$output_dir/gatk.g.vcf.gz" ]; then
    gatk --java-options "-Xmx4g" HaplotypeCaller \
      -L "$peaks" \
      -R "$genome" \
      -I "$bam" \
      -O "$output_dir/gatk.g.vcf.gz" \
      -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation
fi

# correct the genotypes that come out of haplotype caller
if [ ! -f "$output_dir/gatk.vcf.gz" ]; then
    gatk --java-options "-Xmx4g" GenotypeGVCFs \
      -all-sites \
      -L "$peaks" \
      -R "$genome" \
      -V "$output_dir/gatk.g.vcf.gz" \
      -O "$output_dir/gatk.vcf.gz"
fi

gatk --java-options "-Xmx4g" SelectVariants \
  -L "$peaks" \
  -R "$genome" \
  -V "$output_dir/gatk.vcf.gz" \
  --select-type INDEL --select-type NO_VARIATION \
  -O "$output_dir/gatk-snps.vcf.gz"

# convert vcf to table
gatk --java-options "-Xmx4g" VariantsToTable \
  -V "$output_dir/gatk-snps.vcf.gz" \
  -O "$output" \
  -L "$peaks" \
  -F CHROM -F POS -F REF -F ALT -F QD -F FS \
  -F SOR -F MQ -F ReadPosRankSum -GF DP -GF GT;
