#!/usr/bin/env bash

# param1+ (the filtering expression): specify the column name, a comparison operator (>, <, ==), and the value to compare against (ex: 'gatk-indel~DP>20')
# The tab delimited table you want filtered must be passed via stdin and is output via stdout
# Note that the comparison operator and value are fed directly to awk, so make sure to quote them appropriately
# You can specify multiple filtering expressions each as separate parameters to the script

# Currently supported comparison operators are below:
ops=">|<|=="
# But theoretically, you could use any comparison operators supported by awk; just add them to the pipe delimited list above

# TODO: investigate whether we can use miller (http://johnkerl.org/miller) instead of this script

# are there columns to filter on?
if [ ! -z "$1" ]; then
	# read the header into a variable
	read -r head
	echo "$head"
		# transform the provided parameters into awk code:
		# 1) use printf to write each parameter on a separate line
		# 2) use sed to extract just the column names
		# 3) for each colname, determine its numerical index among the other cols
		# 4) use paste to merge the column indices with the text following the ops
		# 5) use grep to ignore columns that aren't in the table
		# 6) use sed to prefix every column idx with a dollar sign $
		# 7) use paste to merge all filtering expressions into a single str
		# 8) use sed to add " && " between each filtering expression
	awk_code="$(
		printf '%s\n' "$@" | sed -r 's/('"$ops"').*$//' | {
			while IFS= read -r colname; do
				printf '%s\n' "$(tr '\t' '\n' <<< "$head" | grep -Fn "$colname" | cut -f1 -d: | head -n1)"
			done
		} | paste -d'\0' - <(printf '%s\n' "$@" | grep -Eo '('"$ops"').*$') | \
		grep -Ev '^('"$ops"')' | sed 's/^/\$/' | paste -s -d $'\t' | sed -r 's/\t/ \&\& /g'
	)"
	# use awk to filter the table
	awk -F $'\t' -v 'OFS=\t' "$awk_code" -
else
	cat
fi
