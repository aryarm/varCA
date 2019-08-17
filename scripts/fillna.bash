#!/usr/bin/env bash

# this script replaces NA values in a large TSV with the provided replacements

# param1: the .tsv.gz file containing the large table for which we want to fill NA values
# param2+: the name of a column followed by the NA value replacement for that col (as two adjacent arguments)
# 		   you can specify multiple column-value pairs
# return (to stdout): the same table except where all specified NA values have been replaced


script_dir="$(dirname "$BASH_SOURCE")";

pfiles=()
for (( i=2; i<$#; i+=2 )); do
	j=$((i+1))
	reg="${!i}"
	val="${!j}"
	# retrieve the idx of the column among the columns
	reg="$(zcat "$1" | head -n1 | tr '\t' '\n' | grep -Fn "$reg" | cut -f1 -d: | head -n1)"
	pfiles+=("{if (\$$reg==\".\") \$$reg=\"$val\"}")
done

zcat "$1" | awk -F $'\t' -v 'OFS=\t' "${pfiles[*]} 1"
