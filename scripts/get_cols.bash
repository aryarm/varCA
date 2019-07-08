#!/bin/bash

# return cols from a file where the column names match a regex filter
# specify the regex filter in arg1 and the file in arg2
# input is read from stdin if arg2 is not specified
# file is assumed to be delimited by tabs
# output is written to stdout

# example to retrieve CHROM, POS, and all ALT cols from a gzipped file:
# ./get_cols.bash '^CHROM$|^POS$|.*~ALT$' <(zcat file.tsv.gz)

col_idxs() {
	# get comma separated list of column indices
	echo "$2" | tr '\t' '\n' | grep -nE "$1" | cut -f1 -d: | tr '\n' ',' | sed 's/,$//'
}

if test -n "$2"; then # work with the file
	head="$(head -n1 "$2")"
	cols="$(col_idxs "$1" "$head")"
	test -z "$cols" && exit 1
	cut -f "$cols" "$2"
elif test ! -t 0; then # work with stdin
	read -r head
	cols="$(col_idxs "$1" "$head")"
	test -z "$cols" && exit 1
	cat <(echo "$head") - | cut -f "$cols"
else
	echo "No data provided..." 2>&1
	exit 1
fi
