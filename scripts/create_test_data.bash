#!/usr/bin/env bash

# Create a small test dataset from a larger example dataset composed of chr1 from the
# Jurkat and MOLT4 cell lines


create_seq() {
    # arg1 (string): either "jurkat" or "molt4"

    num_peaks=2
    if [ "$1" == "molt4" ]; then
        num_peaks=3
    fi

    comm -12 <(
        bedtools intersect -c -bed -a out/peaks/"$1"_chr1/"$1"_chr1_peaks.narrowPeak -b <(
            zcat out/merged_indel/"$1"_chr1/merge.tsv.gz | \
            scripts/cgrep.bash - -E '^(CHROM|POS)$|(gatk|varscan|vardict|pindel|illumina-strelka|pg-indel).*~(ALT)$' | \
            awk -F $'\t' -v 'OFS=\t' '$3 == $4 && $4 == $5 && $5 == $6 && $6 == $7 && $7 == $3 && $3 != "."' | \
            tail -n+2 | \
            awk -F $'\t' -v 'OFS=\t' '{print $1, $2, $2+length($3);}'
        ) | \
        grep -vE '0$' | cut -f-3 | sort
    ) <(
      bedtools intersect -c -bed -a out/peaks/"$1"_chr1/"$1"_chr1_peaks.narrowPeak -b <(
          zcat out/merged_snp/"$1"_chr1/merge.tsv.gz | \
          scripts/cgrep.bash - -E '^(CHROM|POS)$|(gatk|varscan|vardict|pg-snp).*~(ALT)$' | \
          awk -F $'\t' -v 'OFS=\t' '$3 == $4 && $4 == $5 && $3 != "."' | \
          tail -n+2 | \
           awk -F $'\t' -v 'OFS=\t' '{print $1, $2, $2+length($3);}'
        ) | \
        grep -vE '0$' | cut -f-3 | sort
    ) | sort -t $'\t' -k1,1 -k2,2n | head -"$num_peaks" > data/"$1".mini.bed

    samtools sort -n out/align/"$1"_chr1/rmdup.bam | \
    bedtools pairtobed -abam - -b data/"$1".mini.bed 2>/dev/null | \
    samtools view -h | sed 's/'"$1"'_chr1/'"$1"'_mini/g' | grep -Ev '^\@PG' | \
    samtools sort -nO bam > data/"$1"_mini.bam

    if [ "$1" == "molt4" ]; then
        samtools fastq -1 data/"$1".mini.1.fq.gz -2 data/"$1".mini.2.fq.gz -0 /dev/null -s /dev/null -n data/"$1"_mini.bam
        rm data/molt4_mini.bam
    fi
}

create_seq jurkat
create_seq molt4

rm data/hg38.mini.*

cat data/*.bed | sort -k1,1 -k2,2n | sed -e 1b -e '$!d' | tail -n+2 | \
awk -F $'\t' -v 'OFS=\t' '{print $1, 0, $3;}' | \
bedtools slop -i - -g <(cut -f-2 out/data-full-chr1/hg38.chr1.fa.fai) -b 500 | \
bedtools getfasta -fi out/data-full-chr1/hg38.chr1.fa -bed - | \
sed 's/>chr1:0-.*$/>chr1/' > data/hg38.mini.fa

samtools faidx data/hg38.mini.fa
gatk CreateSequenceDictionary -R data/hg38.mini.fa
bwa index data/hg38.mini.fa &>/dev/null

samtools view -h data/jurkat_mini.bam | sed 's/LN:.*$/LN:'"$(cut -f2 data/hg38.mini.fa.fai)"'/' | samtools sort -O bam > data/jurkat.mini.bam
samtools index data/jurkat.mini.bam

rm data/molt4.mini.bed data/jurkat_mini.bam

gatk ValidateSamFile -R data/hg38.mini.fa -I data/jurkat.mini.bam --MODE SUMMARY

