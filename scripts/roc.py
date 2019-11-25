#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import auc
from itertools import product


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
    parser.error("Specify the path to at least one space separated table (w/o a header) with two rows (recall/precision) using options like --gatk-indel path/to/gatk-table")
args = parser.parse_args()
all_args = vars(args)


def get_marker():
    """ retrieve a unique marker in a deterministic order """
    yield from product(range(2, 6), range(1, 3), [52])
    # if that wasn't enough, create some more just at a different angle
    yield from product(range(2, 8), range(1, 3), [0])
    # if we ever get to this point, we need more than 20 markers
    yield from product([2]+list(range(3, 6, 2)), range(1, 3), [90])


colors = {}
markers = get_marker()
# go through each table and get its name from all of the args
for arg in sorted(all_args.keys()):
    if arg not in known_args:
        table = np.loadtxt(all_args[arg])
        # fpr: 1st row, tpr: 2nd row
        if table.ndim != 1:
            # area = auc(table[0], table[1])
            # colors[arg] = plt.plot(
            #     table[0], table[1],
            #     label=arg+": area={0:0.2f}".format(area)
            # )[0].get_color()
            colors[arg] = plt.plot(
                table[0], table[1],
                label=arg.replace('breakca', 'varca')
            )[0].get_color()
        else:
            area = table[1]
            # check if this pt has a curve of the same name
            # to make sure they're the same color
            if arg.endswith("_pt") and arg[:-3] in colors:
                extra = {'color': colors[arg[:-3]]}
            else:
                extra = {}
            plt.plot(
                table[0], table[1], marker=next(markers), ms=12,
                label=arg.replace('breakca', 'varca')+": height={0:0.2f}".format(area), **extra
            )

plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize='small')
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.ylim([0.0, 1.0])
plt.xlim([0.0, 1.0])
plt.savefig(args.out, bbox_inches='tight', pad_inches=0.5)
