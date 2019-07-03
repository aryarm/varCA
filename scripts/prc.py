#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import auc


parser = argparse.ArgumentParser()
parser.add_argument(
    "out", nargs="?", default=sys.stdout, help="the filename to save the data to"
)
known_args, unknown_args = parser.parse_known_args()
# dynamically parse whatever options the user passes us
count = 0
for arg in unknown_args:
    if arg.startswith('--'):
        parser.add_argument('--{}'.format(arg[2:]), type=argparse.FileType('r'))
        count += 1
if count < 1:
    parser.error("Specify the path to at least one two column (recall/precision) table (w/o a header) using options like --gatk-indel path/to/gatk-pr-file")
args = parser.parse_args()
all_args = vars(args)


# go through each table and get its name from all of the args
for arg in all_args:
    if arg not in known_args:
        table = np.loadtxt(all_args[arg])
        # recall: 1st row, precision: 2nd row
        if table.ndim != 1:
            area = auc(table[0], table[1])
            plt.step(table[0], table[1], where='post', label=arg+": area={0:0.2f}".format(area))
        else:
            area = table[1]
            plt.plot(table[0], table[1], 'o', label=arg+": height={0:0.2f}".format(area))

plt.legend()
plt.xlabel('Recall')
plt.ylabel('Precision')
plt.ylim([0.0, 1.0])
plt.xlim([0.0, 1.0])
plt.title('Precision-Recall')
plt.savefig(args.out, bbox_inches='tight', pad_inches=0.5)
