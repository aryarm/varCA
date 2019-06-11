#!/bin/bash

bam="$1"
[[ -z "$1" ]] && { echo "Parameter 1 is empty" 1>&2; exit 1; }
peaks="$2"
[[ -z "$2" ]] && { echo "Parameter 2 is empty" 1>&2; exit 1; }
genome="$3"
[[ -z "$3" ]] && { echo "Parameter 3 is empty" 1>&2; exit 1; }
output_dir="$4"
[[ -z "$4" ]] && { echo "Parameter 4 is empty" 1>&2; exit 1; }
samp="$5"
[[ -z "$5" ]] && { echo "Parameter 5 is empty" 1>&2; exit 1; }


vardict-java -G "$genome" -N "$samp" -b "$bam" -z -c 1 -S 2 -E 3 "$peaks" | \
# vardict teststrandbias
