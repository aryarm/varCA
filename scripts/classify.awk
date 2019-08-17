#!/usr/bin/env awk -f

# This awk script reads a two column table containing REF and ALT columns and outputs a single column classifying each variant
# It is assumed that the REF and ALT alleles have been normalized (ie trimmed and left-aligned), that breakends have been removed, and that multiallelics have been split
# You may optionally pass the name to use for the generated column as a variable called "colname"
# Note that variants that are of "MIXED" type (containing a mix of indels and SNPs) are currently categorized as MNP
# this algorithm was adapted from https://genome.sph.umich.edu/wiki/Variant_classification#Classification_Procedure

# ex usage: ./classify.awk -v 'colname=TYPE' REF-ALT.tsv
# ex-usage for filtering a TSV: cut -f 3,4 CHROM-POS-REF-ALT.tsv | ./classify.awk | grep -En 'INS|DEL|\.' | sed 's/:.*$/p/' | sed -n -f - CHROM-POS-REF-ALT.tsv


BEGIN {if (length(colname) != 0) print colname}
{f=0}
# match regular variants only
$1~/^[ACGTacgt]*$/ && $2~/^[ACGTacgt]*$/ {
	f=1;
	l=length($2)-length($1);
	if (l == 0) {
		if (length($1) == 1 && $1 != $2) {
			print "SNP"
		} else {
			print "MNP"
		}
	} else {
		if (l > 0 && index($2, $1) == 1) {
			print "INS"
		} else if (l < 0 && index($1, $2) == 1) {
			print "DEL"
		} else {
			print "MNP"
		}
	}
}
# match any non-variants
$1~/^[ACGTacgt\.]*$/ && $2~/^[\.]*$/ {
	f=1;
	print "."
}
# match any SVs
$1~/^[ACGTacgt]*$/ && $2~/^<.*>$/ {
	f=1;
	print substr($2, 2, length($2)-2)
}
# match anything else
f==0{print "."}


# EXAMPLES / TEST CASES
# you can check whether the code passes the test cases with the following command:
# tac classify.awk | awk '{if(/^## /)exit;else print}' | sed 's/^# //' | { test="$(tac)"; diff -ys <(echo "$test" | cut -f 3) <(echo "$test" | cut -f 1,2 | ./classify.awk); }
## REF	ALT	TYPE
# A	G	SNP
# AG	A	DEL
# A	AG	INS
# AG	TC	MNP
# A	<DEL>	DEL
# A	<DUP>	DUP
# A	<INV>	INV
# A	.	.
# .	.	.
# .	A	.
# AC	G	MNP
# A	GC	MNP
# AAG	AA	DEL
# AAG	GA	MNP
# AA	AAG	INS
# AG	GAA	MNP
