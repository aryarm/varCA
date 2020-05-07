#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve


parser = argparse.ArgumentParser()
parser.add_argument(
    "-o", "--out", type=argparse.FileType('w', encoding='UTF-8'),
    default=sys.stdout, help="the filename to save the data to"
)
parser.add_argument(
    "-s", "--sorted", action='store_true', help="whether the data is already sorted by its probabilities; the second column will be ignored"
)
parser.add_argument(
    "-f", "--flip", action='store_true', help="whether to flip the probabilities"
)
parser.add_argument(
    "--flip--sorted", action='store_true', help="whether to only flip the probabilities (this is the same as -f and -s specified together and overrides the others)"
)
parser.add_argument(
    "-r", "--roc", action='store_true', help="create roc (instead of prc) data"
)
parser.add_argument(
    "-t", "--thresh", action='store_true', help="also output thresholds for each precision/recall value"
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
    index_col=False, dtype={'probs': np.float_, 'truth': np.bool_},
    low_memory=False, na_values='.', usecols=['truth', 'probs'][:2-args.sorted]
)
df.fillna(0, inplace=True)


if args.sorted:
    print("Input already sorted.", file=sys.stderr)
    scores = df.index/df.index[-1]
    # predictions are already technically inverted
    if not args.flip:
        scores = 1-scores
    else:
        print("Inverting predictions.", file=sys.stderr)
else:
    # replace inf values with a number 1 larger than the next largest value
    if df['probs'].max() == np.float_('inf'):
        df['probs'] = df['probs'].replace(
            np.float_('inf'), np.sort(df['probs'].unique())[-2]+1
        )
    # turn the scores into probabilities if they're not already
    scores = df['probs']/df['probs'].max()
    if args.flip:
        print("Inverting predictions.", file=sys.stderr)
        scores = 1-scores
if args.roc:
    fpr, tpr, thresh = roc_curve(df['truth'], scores)
    if args.thresh:
        np.savetxt(args.out, np.array([fpr, tpr, np.hstack(([0], thresh))]))
    else:
        np.savetxt(args.out, np.array([fpr, tpr]))
else:
    precision, recall, thresh = precision_recall_curve(df['truth'], scores)
    if args.thresh:
        np.savetxt(args.out, np.array([recall, precision, np.hstack(([0], thresh))]))
    else:
        np.savetxt(args.out, np.array([recall, precision]))
