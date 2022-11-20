#!/usr/bin/env bash

# this script replaces NA values in a large TSV with the provided replacements

# stdin: the .tsv file containing the large table for which we want to fill NA values
# param1+: the name of a column followed by the NA value replacement for that col (as two adjacent arguments)
# 		   you can specify multiple column-value pairs
# return (to stdout): the same table except where all specified NA values have been replaced


read -r head

pfiles=()
for (( i=1; i<$#; i+=2 )); do
	j=$((i+1))
	reg="${!i}"
	val="${!j}"
	# retrieve the idx of the column among the columns
	reg="$(echo "$head" | tr '\t' '\n' | grep -Fn "$reg" | cut -f1 -d: | head -n1)"
	pfiles+=("{if (\$$reg==\".\") \$$reg=\"$val\"}")
done

{ echo "$head"; cat; } | awk -F $'\t' -v 'OFS=\t' "${pfiles[*]} 1"
