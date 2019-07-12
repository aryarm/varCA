#!/bin/sh
# help
if [ "$1" == "-h" ] ; then
    echo "Usage= ./breakCA.bash -a <path to R> -b <.bam> -p <.bed> -o <output directory> -g <fasta file for genome>"
    exit 0
fi

# running with options
while getopts ":a:b:p:o:g:" opt
   do
     case $opt in
		a ) path=$OPTARG;;
        b ) bam=$OPTARG;;
        p ) peaks=$OPTARG;;
		o ) output=$OPTARG;;
		g ) genome=$OPTARG;;
     esac
done

# make output directory
mkdir -p $output

# get reads as a tsv file
$path/Rscript --vanilla bin/get_reads_from_bam.R $bam $peaks $output/read.tsv

# get read pileups; allow for correction for mis-alignments
samtools mpileup -B -f $genome -l $peaks $bam | gzip > $output/read.pileups.gz

# get bases written into a file to use during counting
less $output/read.pileups.gz| awk -v OFS='\t' '{print $1,$2,$5}' > $output/read.pileup

# get insertion containing positions
less $output/read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /\+[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > $output/insertion.pileups

# get deletion containing positions
less $output/read.pileups.gz | awk -v OFS='\t' '{ if ($5 ~ /-[0-9]+[ACGTNacgtn]+/) print $1,$2,$5}' > $output/deletion.pileups

# counts reads at each bp position within peaks
$path/Rscript --vanilla bin/count_reads_per_base.R $output/read.tsv $output/insertion.pileups $output/deletion.pileups $output/read.pileup $output/sc.tsv $output/counts.tsv

# add clipping information including information content and clip length
$path/Rscript --vanilla bin/get_clipping_information.R $output/sc.tsv $output/sc_w_seq.tsv $output/clip.info.txt

# calculate posterior mean and standard deviation for clipped positions
$path/Rscript --vanilla bin/calculate_posterior.R $output/counts.tsv $output/posteriors.tsv $output/clip.info.txt $output/all.positions.tsv

echo "done"
