#!/usr/bin/env python3

import argparse
from sys import stdout
import dask.dataframe as dd
from functools import reduce

parser = argparse.ArgumentParser()
parser.add_argument("-k", "--key", default="CHROM\tPOS",
    help="whitespace separated column names that are common to every file")
# parser.add_argument("-o", "--output", required=True,
#     help="the path of the file to output, containing merged columns")
parser.add_argument("files", nargs="*",
    help="csv files whose columns you'd like merged")
args = parser.parse_args()

# 1) Read in the csv files.
# 2) Merge the csv files.
# 3) Write the output.
reduce(
    lambda df1, df2: df1.merge(df2, how='outer', on=args.key.split()),
    (
        dd.read_csv(
            file,
            sep="\t",
            dtype='str',
            blocksize=None,
            compression='gzip' if file.endswith('.gz') else 'infer'
        )
        for file in args.files
    )
).fillna('.').compute().set_index(args.key.split()).to_csv(
    stdout, index=True, na_rep='.', sep='\t'
)
