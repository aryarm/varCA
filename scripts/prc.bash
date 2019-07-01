#!/bin/bash
# param 1: the large tsv
# param 2: the caller id
# param 3: the column to threshold on or "(REF|ALT)" if the type of each variant should be used (note that using "(REF|ALT)" forces param 8 to be 1)
# param 4: the truth set caller id
# param 5: which type of variant to create plots for (ex: DEL, INS, SNP, .); separate each variant by commas if you'd like to micro average over multiple of them
# param 6: the output file or '-' if stdout
# param 7 (optional): whether to invert the threshold column (1) or leave it be (0) (default is 0)
# param 8 (optional): whether to handle the sorting internally (1) or let python do it (0)
# make sure that the large tsv has columns in the order caller~REF, caller~ALT, caller~sort, truth~REF, truth~ALT


script_dir="$(dirname "$BASH_SOURCE")";

function binarize() {
	# use awk to convert each column to binary 1s or 0s based on what type of
	# variant we want to create plots for

	# if the user specified '.' as the variant, label '.' as 0 and everything else as 1
	if [ "$1" == '.' ]; then
		awk '!/\./ {print 1} /\./ {print 0}'
	# if the user specified a comma-separated list of variant types, aggregate the resuts of binarizing all of them
	elif [[ $1 == *","* ]]; then
		local col="$(cat)"
		# create an array of variant types called "$variants"
		IFS=',' read -ra variants <<< "$1"
		# call binarize on the same input with each variant type
		for variant in "${variants[@]}"; do
			echo "$col" | binarize "$variant"
		done
	# if the user specified 'INS', 'DEL', or 'SNP' as the variant, label it as 1 and everything else as 0
	else
		awk '/'"$1"'/ {print 1} !/'"$1"'/ {print 0}'
	fi
}

# use paste to create a two column table of true variants (col1) and their predictions (col2)
# if needed, we binarize the columns according to what the user has defined as the positive label
paste <(
	zcat "$1" | "$script_dir"/get_cols.bash "^$4~(REF|ALT)$" | tail -n+2 | \
	"$script_dir"/classify.awk | binarize "$5"
) <(
	zcat "$1" | "$script_dir"/get_cols.bash "^$2~$3$" | tail -n+2 | {
		# if the predictor is specified as "(REF|ALT)", we create a new column
		# containing the variant types and then binarize them, just like the
		# truth column
		# otherwise, we expect that the column contains numerical values
		if [ "$3" == "(REF|ALT)" ]; then
			"$script_dir"/classify.awk | binarize "$5"
		else
			# make sure to repeat the col n+1 times where n is the number of
			# commas in the variant type string "$5"
			col="$(cat)"
			for (( i = 0; i < "$(awk -F, '{print NF}' <<< "$5")"; i++ )); do
				echo "$col"
			done
		fi
	}
) | {
	# if "(REF|ALT)" is specified as the prediction column, the predictions
	# must be sorted internally because sklearn will not perform a 'stable'
	# sort (ie one in which ties are broken by using the original ordering)
	# unstable sorts lead to very rectangular-ish curves
	if [ "$8" == "1" ] || [ "$3" == "(REF|ALT)" ]; then
		sort -s -t $'\t' -k2,2nr
	else
		cat
	fi
} | \
python "$script_dir"/statistics.py -o "$6" "$([ "$7" == "1" ] && echo "--flip")$([ "$8" == "1" ] || [ "$3" == "(REF|ALT)" ] && echo "--sorted")";
