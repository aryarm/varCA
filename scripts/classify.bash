#!/bin/bash

# param1: the .tsv.gz file containing the large table we want to classify
# param2: which type of variant to use as the positive label (ex: DEL, INS, SNP, .); separate each variant by commas if you'd like to output multiple labels or pass the empty string '' if you want all labels
# return (to stdout): the same table, where all REF and ALT columns have been transformed into categorical classification columns (not in the original order)


script_dir="$(dirname "$BASH_SOURCE")";

# get an array of the ref and alt column indices for each caller
# ex: if gatk-indel's REF and ALT columns appear in columns 4 and 5 respectively, the following would be an entry in the array: "gatk-indel:4,5"
callers=($(
	zcat "$1" | head -n1 | tr '\t' '\n' | \
	grep -Pno '^.*~(REF|ALT)$' | tr : '\t' | tr '~' '\t' | \
	sort -k2,2 -k3,3r | cut -f -2 | \
	awk -F '\t' -v 'ORS= ' 'NR%2{printf "%s:%s,",$2,$1;next;} {print $1}'
))

# retrieve a comma separated string of all the REF/ALT col indices
all_idx="$(echo "${callers[@]}" | tr ' ' '\n' | sed 's/.*://' | paste -s -d,)"

function binarize() {
	# use awk to convert each column to binary 1s or 0s based on what type of
	# variant we want to use as the pos label
	
	# first, rename the header
	read -r head

	# if the user didn't specify a type of variant, don't do any filtering
	if [ -z "$1" ]; then
		echo "$head:$1"
		cat
	# if the user specified '.' as the variant, label '.' as 0 and everything else as 1
	elif [ "$1" == '.' ]; then
		echo "$head:$1"
		awk '!/\./ {print 1} /\./ {print 0}'
	# if the user specified a comma-separated list of variant types, aggregate the resuts of binarizing all of them
	elif [[ $1 == *","* ]]; then
		local col="$(cat)"
		# create an array of variant types called "$variants"
		IFS=',' read -ra variants <<< "$1"
		# call binarize on the same input with each variant type
		for i in "${!variants[@]}"; do
			variants[$i]="<(cat <(echo \"\$head\") <(echo \"\$col\") | binarize "${variants[$i]}")"
		done
		eval "paste ${variants[@]}"
	# if the user specified 'INS', 'DEL', or 'SNP' as the variant, label it as 1 and everything else as 0
	else
		echo "$head:$1"
		awk '/'"$1"'/ {print 1} !/'"$1"'/ {print 0}'
	fi
}

# iterate through the REF,ALT column indices of each caller
for i in "${!callers[@]}"; do
	id="$(cut -f1 -d: <<< "${callers[$i]}")"
	idx="$(cut -f2 -d: <<< "${callers[$i]}")"
	callers[$i]="<(zcat \"\$1\" | cut -f "$idx" | \"\$script_dir\"/classify.awk -v \"colname="$id"~CLASS\" | binarize \"\$2\")"
done

eval "paste ${callers[@]}"
