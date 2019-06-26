#!/usr/bin/env python3

import sys
# import gzip
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

df = pd.read_csv(
    args.table, sep=',', header=None, names=['caller', 'truth']
    index_col=False, dtype=str, low_memory=False, na_values='.'
)


categories = df['caller'].unique()
rates = {
    category: {
        'TP': [],
        'FP': [],
        'TN': [],
        'FN': []
    } for category in categories
}


thresh = 0
curr_class = None
for line in df:
    if curr_class is not None and curr_class != line['truth']:
        # we've hit a threshold
        # output rates
        thresh += 1
    for category in categories:
        if category == line['truth'] and category == line['caller']:
            ctype = 'TP'
        elif category == line['truth']:
            ctype = 'FP'
        elif category == line['caller']:
            ctype = 'FN'
        else:
            ctype = 'TN'
        # check if a new treshold number exists before adding to it
        if thresh == rates[category][ctype].size:
            rates[category][ctype].append(0)
        # add to the current ctype number
        rates[category][ctype][thresh] += 1
    curr_class = line['truth']

totals = {
    category: {
        ctype: rates[category][ctype][-1]
        for ctype in rates[category]
    } for category in categories
}

# calculate precision recall by micro-averaging for each class (at each threshold)
for i in range(thresh+1):
    TP_sum = 0
    FP_sum = 0
    FN_sum = 0
    for category in categories:
        TP_sum += totals[category]['TP'][thresh]
        FP_sum += totals[category]['FP'][thresh]
        FN_sum += totals[category]['FN'][thresh]
    rates['precision'].append(TP_sum/(TP_sum+FP_sum))
    rates['recall'].append(TP_sum/(TP_sum+FN_sum))

plt.scatter(rates['recall'], rates['precision'])
plt.show()
