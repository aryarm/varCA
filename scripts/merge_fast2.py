#!/usr/bin/env python3

from sys import stdout
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("-b", "--bed", required=True,
    help="the path of the bed file containing the peaks")
parser.add_argument("files", nargs="*",
    help="csv files whose columns you'd like merged")
args = parser.parse_args()

# parse the bed file
for peak in args.bed:
    peak.split('\t')

# output.close()
# for file in files:
#     file.close()
