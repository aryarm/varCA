#!/bin/bash

# param1: the .tsv.gz file containing the large table we want to classify
# param2+: a regex pattern for a column followed by the NA value replacement for those cols (separated by a space)
# 		   you can specify multiple column-value pairs like this
# return (to stdout): the same table, where all REF and ALT columns have been transformed into categorical classification columns (not in the original order)


script_dir="$(dirname "$BASH_SOURCE")";

pfiles=()
for (( i = 2; i < $#; i+=2 )); do
	j=$((i+1))
	reg="${!i}"
	val="${!j}"
	pfiles+=("<(zcat \"\$1\" | \"\$script_dir\"/get_cols.bash '$reg' | awk -F '\t' -v 'OFS=\t' '{{for (i=1; i<=NF; i++) if (\$i==\".\") \$i=\"$val\"}}1')")
done

eval "paste ${pfiles[@]}"
