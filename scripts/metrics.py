#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
import sklearn.metrics


parser = argparse.ArgumentParser()
parser.add_argument(
    "-o", "--out", default=sys.stdout, help="the filename to save the data to"
)
parser.add_argument(
    "-m", "--metrics", default='p,r,f', help="a comma separated list of metrics to output; use 'p' for precision, 'r' for recall, and 'f' for the F-beta score"
)
parser.add_argument(
    "table", nargs="?", default=sys.stdin,
    help="a two column (truth/predicted) table of variant classifications w/o a header"
)
args = parser.parse_args()

# read the file into a pandas data frame
df = pd.read_csv(
    args.table, sep='\t', header=None, names=['truth', 'predict'],
    index_col=False, dtype=str,
    low_memory=False, na_values='.'
)

to_idx = {'p': 0, 'r': 1, 'f': 2}
metrics = [to_idx[metric] for metric in args.metrics.split(",")]
scores = np.array(
    sklearn.metrics.precision_recall_fscore_support(df['truth'], df['predict'])
)

np.savetxt(args.out, scores[metrics])
