#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import pandas as pd
import sklearn.metrics


parser = argparse.ArgumentParser()
parser.add_argument(
    "-o", "--out", type=argparse.FileType('w', encoding='UTF-8'),
    default=sys.stdout, help="the filename to save the data to"
)
parser.add_argument(
    "-m", "--metrics", default='r,p,b,t,f', help=(
        "a comma separated, ordered list of metrics to output; use 'p' for "
        "precision, 'r' for recall, 'b' for the F-beta score, 't' for total "
        "positives, and 'f' for total negatives"
    )
)
# parser.add_argument(
#     "-n", "--names", action='store_true', help="whether to include the names of the metrics alongside them"
# )
parser.add_argument(
    "table", nargs="?", default=sys.stdin,
    help="a two column (truth/predicted) table of variant classifications w/o a header"
)
args = parser.parse_args()

# read the file into a pandas data frame
df = pd.read_csv(
    args.table, sep='\t', header=None, names=['truth', 'predict'],
    index_col=False, dtype={'predict': np.bool_, 'truth': np.bool_},
    low_memory=False, na_values='.'
)

to_idx = {'p': 0, 'r': 1, 'b': 2, 'f': 4, 't': 5}
metrics = [to_idx[metric] for metric in args.metrics.split(",")]
scores = np.append(
    sklearn.metrics.precision_recall_fscore_support(
        df['truth'], df['predict'], average='binary'
    ), sklearn.metrics.confusion_matrix(
        df['truth'], df['predict']
    ).sum(0)
)

np.savetxt(args.out, scores[metrics], fmt='%f')
