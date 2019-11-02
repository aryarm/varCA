#!/usr/bin/env awk -f

# This awk script attempts to normalize numerical values in tabular data. Run tail on this file for some examples.

BEGIN{
	CONVFMT="%.10f"
}
{
	for (i=1; i<=NF; i++) {
		if ($i ~ /[+-]?(?:0|[1-9][0-9]*)(?:\.[0-9]*)?[eE][+-]?[0-9]+/) {
			$i+=0
		}
		if ($i ~ /^[0-9]+(\.[0-9]+)?%$/) {
			$i=substr($i, 1, length($i)-1)
		}
	}
}
1

# EXAMPLES / TEST CASES
# you can check whether the code passes the test cases with the following command:
# tac norm_nums.awk | awk '{if(/^## /)exit;else print}' | sed 's/^# //' | { test="$(tac)"; diff -ys <(echo "$test" | cut -f 2) <(echo "$test" | cut -f 1 | ./norm_nums.awk); }
## NUM	RESULT
# 0.23	0.23
# a str	a str
# 4%	4
# -4.8e-2%	-0.0480000000
# 2.5e5%	250000
# +2.500004e+5%	250000.4000000000
