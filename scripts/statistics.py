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


# # extract the possible classes
# # (ex: INS, DEL, NA for a dataset of indels)
# classes = df['caller'].unique()
# # what are the test statistics?
# # true positives (TP), false positives (FP), false negatives (FN), and true
# # negatives (TN)
# statistics = ['TP', 'FP', 'FN', 'TN']

# # initialize a data structure for storing the test statistics at each threshold
# rates = np.zeros(
#     (len(classes), len(df.index), len(statistics)),
#     dtype='uint32'
# )


# # start calculating the test statistics
# thresh = 0
# curr_class = None
# for idx, line in df.iterrows():
#     if curr_class is not None and curr_class != tuple(line):
#         # we've hit a threshold
#         thresh += 1
#         # make sure to copy the statistics from the previous row for each class
#         for category in range(len(classes)):
#             rates[category][thresh] = rates[category][thresh-1]
#     for category in range(len(classes)):
#         # use binary to decimal conversion to figure out which statistic to +1
#         statistic = int(classes[category] == line['truth']) * 2 + \
#             int(classes[category] == line['caller'])
#         # add to the current statistic number
#         rates[category][thresh][statistic] += 1
#     curr_class = tuple(line)


# # only keep np rows that we used
# rates = rates[:, :thresh+1, :]

# # calculate precision/recall if desired
# # if args.type == 'prc':
# #     # first, precision
# #     classes.append('precision')
# #     rates[:, :thresh+1, len(classes)-1] = np.array()
# #     # then, recall
# #     classes.append('recall')
# #     rates[:, :thresh+1, len(classes)-1] = np.array()

# # write data to stdout
# print("## shape:{0[0]},{0[1]},{0[2]}".format(rates.shape))
# print("## classes:" + ",".join(classes))
# print("")
# print("# total")
# np.savetxt(sys.stdout, rates[:, -1, :])
# print("")
# for category in range(len(classes)):
#     print("# {0}".format(classes[category]))
#     np.savetxt(sys.stdout, rates[category])
#     print("")
