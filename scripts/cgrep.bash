#!/usr/bin/env bash

# return cols from a file where the column names match a regex PATTERN
# usage: ./cgrep.bash FILE [...] PATTERN (where FILE can be - for stdin)
# FILE is assumed to be delimited by tabs
# extra arguments to grep can be provided in [...]
# if stdin is provided without the -, the col idxs will be output instead

# example to retrieve CHROM, POS, and all ALT cols from a gzipped file:
# zcat file.tsv.gz | ./cgrep.bash '^CHROM$|^POS$|.*~ALT$' - | less

col_idxs() {
	# get comma separated list of column indices
	tr '\t' '\n' | grep -n "$@" | cut -f1 -d: | paste -s -d',' | sed 's/,$//'
}

# first, check: should we work with the file or stdin?
if [ -f "$1" ]; then # work with the file
	cols="$(head -n1 "$1" | { shift && col_idxs "$@"; })"
	test -z "$cols" && exit 1
	cut -f "$cols" "$1"
elif [ ! -t 0 ]; then # work with stdin
	read -r head
	cols="$(echo "$head" | { shift && col_idxs "$@"; })"
	test -z "$cols" && exit 1
	if [ "$1" = "-" ]; then
		{ echo "$head"; cat -; } | cut -f "$cols"
	else
		echo "$cols"
	fi
else
	echo "No data provided..." 2>&1; exit 1
fi
