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
samp="$6"
[[ -z "$6" ]] && { echo "Parameter 6 is empty" ; exit 1; }

pindel -f "$genome" -i <(echo "$bam" 300 "$samp") -j "$peaks" -o "$output_dir/"
pindel2vcf -p 