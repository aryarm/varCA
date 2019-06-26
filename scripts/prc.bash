#!/bin/bash

# param 1: the large tsv
# param 2: the caller id
# param 3: the column to threshold on
# param 4: the truth set caller id
# param 5 (optional): whether to sort the column in reverse order (default is true)


# make sure that the large tsv has columns in the order caller~REF, caller~ALT, caller~sort, truth~REF, truth~ALT


zcat "$1" | \
"$(dirname "$BASH_SOURCE")"/get_cols.bash "^$2~(REF|ALT|$3)$|^$4~(REF|ALT)$" | {
	read -r head;
	echo "$head";
	sort -t $'\t' -k3,3n"$(test -z "$5" && echo "r")";
} | {
 	variants="$(cat)";
 	paste -d, \
 	<(
 		echo "$variants" | \
 		"$(dirname "$BASH_SOURCE")"/get_cols.bash "^$2~(REF|ALT)" | \
 		tail -n+2 | "$(dirname "$BASH_SOURCE")"/classify.awk
 	) \
 	<(
 		echo "$variants" | \
 		"$(dirname "$BASH_SOURCE")"/get_cols.bash "^$4~(REF|ALT)" | \
 		tail -n+2 | "$(dirname "$BASH_SOURCE")"/classify.awk
	) | \
 	python "$(dirname "$BASH_SOURCE")"/prc.py
}
