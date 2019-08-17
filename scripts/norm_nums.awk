#!/usr/bin/env awk -f

# This awk script attempts to normalize numerical values in tabular data. Here are some examples:
# 1.75e-1  ==>  0.1750000000
# 4%       ==>  4
# 4.8%     ==>  4.8
# 4.8e-2%  ==>  .0048

BEGIN{
	CONVFMT="%.10f"
}
{
	for (i=1; i<=NF; i++) {
		if ($i ~ /[+-]?(?:0|[1-9]\d*)(?:\.\d*)?[eE][+-]?[0-9]+/) {
			$i+=0
		}
		if ($i ~ /^[0-9]+(\.[0-9]+)?%$/) {
			$i=substr($i, 1, length($i)-1)
		}
	}
}
1
