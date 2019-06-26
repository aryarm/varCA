#!/usr/bin/env python3

import sys
import math
import pandas
import argparse
import numpy as np
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser()
parser.add_argument(
    "table", nargs="?", type=argparse.FileType('r'), default=sys.stdin,
    help="a two column (caller/truth) table of variant classifications"
)
args = parser.parse_args()


# read the file into a pandas data frame
df = pd.read_csv(
    args.table, sep=',', header=None, names=['caller', 'truth'],
    index_col=False, dtype=str, low_memory=False
)


# extract the possible classes
# (ex: INS, DEL, NA for a dataset of indels)
classes = df['caller'].unique()
# what are the test statistics?
# true positives (TP), false positives (FP), false negatives (FN), and true
# negatives (TN)
statistics = ['TP', 'FP', 'FN', 'TN']

# initialize a data structure for storing the test statistics at each threshold
rates = np.zeros(
    (len(classes), len(df.index), len(statistics)),
    dtype='uint32'
)

thresh = 0
curr_class = None
for idx, line in df.iterrows():
    if curr_class is not None and curr_class != tuple(line):
        # we've hit a threshold
        thresh += 1
        # make sure to copy the statistics from the previous row for each class
        for category in range(len(classes)):
            rates[category][thresh] = rates[category][thresh-1]
    for category in range(len(classes)):
        # use binary to decimal conversion to figure out which statistic to +1
        statistic = int(classes[category] == line['truth']) * 2 + \
            int(classes[category] == line['caller'])
        # add to the current statistic number
        rates[category][thresh][statistic] += 1
    curr_class = tuple(line)

# clear some memory
del df
# make sure to delete unused np rows now that we know the number of thresholds
# hopefully python garbage collection will free this memory?
rates = rates[:, :thresh+1, :]

# do other stuff...

# calculate precision recall by micro-averaging for each class (at each threshold)
for i in range(thresh+1):
    TP_sum = 0
    FP_sum = 0
    FN_sum = 0
    for category in classes:
        TP_sum += totals[category]['TP'][thresh]
        FP_sum += totals[category]['FP'][thresh]
        FN_sum += totals[category]['FN'][thresh]
    rates['precision'].append(TP_sum/(TP_sum+FP_sum))
    rates['recall'].append(TP_sum/(TP_sum+FN_sum))

plt.scatter(rates['recall'], rates['precision'])
plt.show()
