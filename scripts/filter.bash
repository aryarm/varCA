#!/bin/bash
# param 1 = the filtering expressions: specify the caller id, a tilde "~", the name of the column to use for filtering, a comparison operator (>, <, ==), and the value to compare against (ex: 'gatk-indel~DP>20'); you can specify multiple filtering expressions by separating them with tabs but the order in which they are given must match that of the tsv

filter_cols() { cat; }
# are there columns to filter on?
if [ ! -z "$1" ]; then
	# retrieve an array of the columns to filter on, sorted according to their order in the tsv
	filter_col="$(tr '\t' '\n' <<< "$1" | sed -r 's/(>|<|==).*$//')"
	# prepare a pipe delimited list of the columns to filter on
	filter_cols="$(paste -s -d'|' <<< "$filter_col")"
	# prepare a function for filtering the rows
	filter_cols() {
		awk -F $"\t" -v 'OFS=\t' "$(gawk -v'RS=,' -v'ORS= && ' '{split($0,a,"(>|<|==)",seps); print "$" ++i+2 seps[length(seps)] a[length(a)] }' <<< "$1" | head -n1)" | cut -f -2;
	}
fi

filter_cols "$1"
