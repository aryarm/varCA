#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve


parser = argparse.ArgumentParser()
parser.add_argument(
    "-o", "--out", type=argparse.FileType('w', encoding='UTF-8'),
    default=sys.stdout, help="the filename to save the data to"
)
parser.add_argument(
    "-s", "--sorted", action='store_true', help="whether the data is already sorted by its probabilities"
)
parser.add_argument(
    "-f", "--flip", action='store_true', help="whether to flip the probabilities"
)
parser.add_argument(
    "--flip--sorted", action='store_true', help="whether to only flip the probabilities (this is the same as -f and -s specified together and overrides the others)"
)
parser.add_argument(
    "table", nargs="?", default=sys.stdin,
    help="a two column (truth/probs) table of variant classifications w/o a header"
)
args = parser.parse_args()
if args.flip__sorted:
    args.flip = True
    args.sorted = True
if args.table == '':
    args.table = sys.stdin

# read the file into a pandas data frame
df = pd.read_csv(
    args.table, sep='\t', header=None, names=['truth', 'probs'],
    index_col=False, dtype={'probs': np.float32, 'truth': np.bool_},
    low_memory=False, na_values='.'
)
df.fillna(0, inplace=True)


if args.sorted:
    print("Input already sorted.", file=sys.stderr)
    scores = df.index/df.index[-1]
    # predictions are already technically already inverted
    if not args.flip:
        scores = 1-scores
    else:
        print("Inverting predictions.", file=sys.stderr)
else:
    scores = df['probs']/df['probs'].max()
    if args.flip:
        print("Inverting predictions.", file=sys.stderr)
        scores = 1-scores
precision, recall, thresh = precision_recall_curve(df['truth'], scores)

np.savetxt(args.out, np.array([recall, precision]))
