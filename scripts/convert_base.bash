#!/bin/bash

###
# convert BED and TSV files between 0 and 1 based coordinate systems
###

# first, check if it's a tsv or bed file
if [[ "$1" =~ *tsv ]]; then
	echo "currently unsupported filetype" 1>&2; exit 1;
elif [[ "$1" =~ *bed ]]; then		
	if [ "$2" == "0" ]; then
		if [ -f "$1" ]; then
			cat "$1" | \
			awk '{printf("%s\t%d\t%d\n",$1,int($2)-1,int($3));}'
		else
			awk '{printf("%s\t%d\t%d\n",$1,int($2)-1,int($3));}'
		fi
	else
		if [ -f "$1" ]; then
			cat "$1" | \
			awk '{printf("%s\t%d\t%d\n",$1,int($2)+1,int($3));}'
		else
			awk '{printf("%s\t%d\t%d\n",$1,int($2)+1,int($3));}'
		fi
	fi
else
	echo "invalid filetype" 1>&2; exit 1;
fi
