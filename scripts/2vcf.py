#!/usr/bin/env python3

import sys
import argparse
import pandas as pd


parser = argparse.ArgumentParser()
parser.add_argument(
    "-o", "--out", default=sys.stdout, help="the filename to save the vcf to"
)
parser.add_argument(
    "caller-vcfs", help="the path to the 'callers' directory for this sample, which should have been craeted by the prepare pipeline"
)
parser.add_argument(
    "results", nargs="?", default=sys.stdin,
    help="a results.tsv.gz file from the output of the classify pipeline"
)
args = parser.parse_args()


