#!/usr/bin/env python3

import sys
import argparse
import matplotlib
import numpy as np
import pandas as pd
from pathlib import Path


parser = argparse.ArgumentParser("Creates a summary metrics.tsv table containing performance metrics for varCA and the variant callers in the ensemble.")
parser.add_argument(
    "pts", type=Path, help="path to the pts dir within your test set's classify/prc/ directory"
)
parser.add_argument(
    "out", nargs='?', default=sys.stdout,
    help="the filename to which to save the metrics table (default: stdout)"
)
args = parser.parse_args()

df = pd.concat([
    pd.read_csv(f, sep="\t", header=None, names=['varca' if f.stem == 'breakca' else f.stem])
    for f in args.pts.rglob('*.txt')
], axis=1)
df.index = ['recall', 'precision', 'f-beta', 'total positives', 'total negatives', 'auroc', 'avg precision']
df.sort_values('f-beta', axis=1, inplace=True, ascending=False)

df.to_csv(args.out, sep="\t")
