#!/usr/bin/env bash

# List all positions (to stdout) at which VarCA called a variant but there were conflicts for the alternative allele.
# Note that this script will only work for the chosen subset of variant callers.

# arg1: merge.tsv.gz file
# arg2: results.tsv.gz file
# arg3: 'snp' or 'indel'

# Ex: scripts/allele_conflicts.bash out/merged_indel/SRR891269/merge.tsv.gz out/new-classify/classify-indel/SRR891269_even_test/results.tsv.gz indel

if [ "$3" = 'indel' ]; then
    zcat "$1" | \
    scripts/cgrep.bash - -E '^(CHROM|POS)$|(gatk|varscan|vardict|pindel|illumina-strelka|pg-indel).*~(ALT)$' | \
    awk -F $'\t' -v 'OFS=\t' '$3 != $4 || $4 != $5 || $5 != $6 || $6 != $7 || $7 != $3' | \
    tail -n+2 | \
    awk -F $'\t' -v 'OFS=\t' '{ { for(i=3; i<=NF-1; i++) count[$i]++; } PROCINFO["sorted_in"] = "@val_num_desc"; { for (val in count) { print $0, val, count[val]; break; } } delete count; }' | \
    sed 's/\t/,/' | \
    LC_ALL=C sort -t $'\t' -k1,1 | \
    LC_ALL=C join -t $'\t' -e '' -j1 -o auto --nocheck-order - <(
        zcat "$2" | \
        awk -F $"\t" -v 'OFS=\t' 'NR == 1 || $NF == 1' | \
        sed 's/\t/,/' | \
        LC_ALL=C sort -t $'\t' -k1,1
    )
else
    zcat "$1" | \
    scripts/cgrep.bash - -E '^(CHROM|POS)$|(gatk|varscan|vardict|pg-snp).*~(ALT)$' | \
    awk -F $'\t' -v 'OFS=\t' '$3 != $4 || $4 != $5 || $5 != $3' | \
    tail -n+2 | \
    awk -F $'\t' -v 'OFS=\t' '{ { for(i=3; i<=NF-1; i++) count[$i]++; } PROCINFO["sorted_in"] = "@val_num_desc"; { for (val in count) { print $0, val, count[val]; break; } } delete count; }' | \
    sed 's/\t/,/' | \
    LC_ALL=C sort -t $'\t' -k1,1 | \
    LC_ALL=C join -t $'\t' -e '' -j1 -o auto --nocheck-order - <(
        zcat "$2" | \
        awk -F $"\t" -v 'OFS=\t' 'NR == 1 || $NF == 1' | \
        sed 's/\t/,/' | \
        LC_ALL=C sort -t $'\t' -k1,1
    )
fi | \
sed 's/,/\t/' | \
LC_ALL=C sort -t $'\t' -k1,1V -k2,2n

# To make a table containing alternative alleles for the following rules (as columns)
# 1) variant caller priority
# 2) majority rule
# 3) platinum genomes
# bcftools query -f '%CHROM\t%POS\t%INFO/CALLER\n' out/new-classify/classify-indel/SRR891269_even_test/final.vcf.gz | sed 's/gatk-indel/1/; s/varscan-indel/2/; s/vardict-indel/3/; s/pindel/4/; s/illumina-strelka/5/' | sed 's/\t/,/' | LC_ALL=C sort -t $'\t' -k1,1 | LC_ALL=C join -t $'\t' -e '' -j1 -o auto --nocheck-order <(zcat out/new-classify/classify-indel/SRR891269_even_test/allele_conflicts.tsv.gz | sed 's/\t/,/' | sort -t $'\t' -k1,1) - | cut -f2- | awk -F $'\t' -v 'OFS=\t' '{print $$NF, $7, $6;}' | less
# bcftools query -f '%CHROM\t%POS\t%INFO/CALLER\n' out/new-classify/classify-snp/SRR891269_even_test/final.vcf.gz | sed 's/gatk-snp/1/; s/varscan-snp/2/; s/vardict-snp/3/' | sed 's/\t/,/' | LC_ALL=C sort -t $'\t' -k1,1 | LC_ALL=C join -t $'\t' -e '' -j1 -o auto --nocheck-order <(zcat out/new-classify/classify-snp/SRR891269_even_test/allele_conflicts.tsv.gz | sed 's/\t/,/' | sort -t $'\t' -k1,1) - | cut -f2- | awk -F $'\t' -v 'OFS=\t' '{print $$NF, $5, $4;}' | less
