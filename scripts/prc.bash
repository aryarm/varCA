#!/bin/bash

# param 1: the large tsv
# param 2: the caller id
# param 3: the column to threshold on
# param 4: the truth set caller id
# param 5 (optional): whether to sort the column in reverse order (default is true)


# make sure that the large tsv has columns in the order caller~REF, caller~ALT, caller~sort, truth~REF, truth~ALT


script_dir="$(dirname "$BASH_SOURCE")"
paste -d $'\t' <(
	zcat "$1" | "$script_dir"/get_cols.bash "^$2~(REF|ALT)$" | \
	tail -n+2 | "$script_dir"/classify.awk
) <(
	zcat "$1" | "$script_dir"/get_cols.bash "^$2~$3$" | \
	tail -n+2
) <(
	zcat "$1" | "$script_dir"/get_cols.bash "^$4~(REF|ALT)$" | \
	tail -n+2 | "$script_dir"/classify.awk
) | \
sort -t $'\t' -k2,2n"$(test -z "$5" && echo "r")" | \
cut -f1,3 | \
python "$script_dir"/prc.py
