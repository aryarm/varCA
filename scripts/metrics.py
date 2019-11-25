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
    "-m", "--metrics", default='r,p,b,t,f,a,v', help=(
        "a comma separated, ordered list of metrics to output; use 'r' for "
        "recall, 'p' for precision, 'b' for the F-beta score, 't' for total "
        "positives, 'f' for total negatives, 'a' for the AUROC, and 'v' for "
        "the avg precision score"
    )
)
# parser.add_argument(
#     "-n", "--names", action='store_true', help="whether to include the names of the metrics alongside them"
# )
parser.add_argument(
    "-p", "--ignore-probs", action='store_true', help="whether to only read truth and predict columns and ignore a probs column if it is provided; note that the AUROC and avg precision will not be output"
)
parser.add_argument(
    "-f", "--flip", action='store_true', help="whether to flip the probabilities; only relevant if --ignore-probs is not passed"
)
parser.add_argument(
    "table", nargs="?", default=sys.stdin,
    help="a three column (truth/predicted/probs) table of variant classifications w/o a header"
)
args = parser.parse_args()

# which cols should we read?
fields = ['truth', 'predict', 'probs']
dtypes = {'predict': np.bool_, 'truth': np.bool_, 'probs': np.float_}
if args.ignore_probs:
    fields = fields[:2]
    dtypes.pop('probs')
# read the file into a pandas data frame
df = pd.read_csv(
    args.table, sep='\t', header=None, names=fields,
    index_col=False, dtype=dtypes,
    low_memory=False, na_values='.'
)
df.fillna(0, inplace=True)

# calculate the metrics
scores = np.append(
    sklearn.metrics.precision_recall_fscore_support(
        df['truth'], df['predict'], beta=0.5, average='binary'
    ),
    sklearn.metrics.confusion_matrix(df['truth'], df['predict']).sum(0)
)
# calculate additional metrics if we can
if not args.ignore_probs:
    # replace inf values with a number 1 larger than the next largest value
    if df['probs'].max() == np.float_('inf'):
        df['probs'] = df['probs'].replace(
            np.float_('inf'), np.sort(df['probs'].unique())[-2]+1
        )
    # turn the scores into probabilities if they're not already
    probs = df['probs']/df['probs'].max()
    if args.flip:
        print("Inverting predictions.", file=sys.stderr)
        probs = 1-probs
    scores = np.append(
        scores,
        np.array([
            sklearn.metrics.roc_auc_score(df['truth'], probs),
            sklearn.metrics.average_precision_score(df['truth'], probs)
        ])
    )

# which metrics should we return?
to_idx = {'p': 0, 'r': 1, 'b': 2, 'f': 4, 't': 5, 'a': 6, 'v': 7}
metrics = [
    to_idx[metric]
    for metric in args.metrics.split(",")
    if not args.ignore_probs or metric not in {'a', 'v'}
]

np.savetxt(args.out, scores[metrics], fmt='%f')
