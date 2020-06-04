#!/usr/bin/env python3
import sys
import pysam
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, recall_score
from sklearn.metrics import roc_curve


# strategy: iterate through each line in the results file and seek to the corresponding line in the merged file
#results = pd.read_csv(args.results, sep="\t", low_memory=True, dtype="str", na_values=['.'], chunksize=100)


# def plot_precision_probs(a):
#   a = np.loadtxt(a)[1:,1:]
#   a[1] = -10*np.log10(1-a[1])
#   p = np.poly1d(np.polyfit(a[0], a[1], 1))
#   plt.scatter(a[0], a[1], label="_nolegend_")
#   plt.xlabel("Precision")
#   plt.ylabel("QUAL")
#   plt.plot(a[0], p(a[0]), label=str(p))
#   plt.legend()

def plot_recall_probs(a):
    df = pd.read_csv(results, sep="\t", header=0, index_col=["CHROM", "POS"], usecols=['CHROM', 'POS', 'breakca~truth', 'breakca~prob.1', 'breakca~CLASS:']).sort_values(by='breakca~prob.1')
    a = df[['breakca~truth', 'breakca~prob.1']].to_numpy().T
    roc = np.array(roc_curve(a[0], a[1]))[:,:-2]
    roc[1] = -10*np.log10(1-roc[1])
    p = np.poly1d(np.polyfit(roc[0], roc[1], 1))
    plt.scatter(roc[0], roc[1], label="_nolegend_")
    plt.xlabel("Recall")
    plt.ylabel("QUAL")
    plt.plot(roc[0], p(roc[0]), label=str(p))
    plt.legend()

def plot_accuracy_probs(results, bin_size):
    """ bin the sites and calculate an accuracy for that bin """
    df = pd.read_csv(results, sep="\t", header=0, index_col=["CHROM", "POS"], usecols=['CHROM', 'POS', 'breakca~truth', 'breakca~prob.1', 'breakca~CLASS:']).sort_values(by='breakca~prob.1')
    lst = np.array([
        (grp.iloc[:,1].mean(), accuracy_score(grp.iloc[:,0], grp.iloc[:,2]))
        for grp in (df.iloc[i:i+bin_size] for i in range(0,len(df)-bin_size+1,bin_size))
    ])
    plt.plot(*lst.T)
    return lst

def phred(val):
    return -10*np.log10(1-val)

def plot_tpr_probs(df, bin_size):
    """ bin the sites and calculate an accuracy for that bin """
    df = df.sort_values(by='breakca~prob.1')
    # split the dataset into high (above 0.5) vs low (below 0.5) scores
    df1 = df[df['breakca~truth']>=0.5]
    df2 = df[df['breakca~truth']<0.5][::-1]
    lst = np.array([
        (grp.iloc[:,1].mean(), phred(recall_score(grp.iloc[:,0], grp.iloc[:,2])))
        for grp in (df1.iloc[i:i+bin_size] for i in range(0,len(df1)-bin_size+1,bin_size))
    ])
    bin_size = int(bin_size*len(df2)/len(df1))
    lst2 = np.array([
        (grp.iloc[:,1].mean(), phred(recall_score(grp.iloc[:,0], grp.iloc[:,2], pos_label=0)))
        for grp in (df2.iloc[i:i+bin_size] for i in range(0,len(df2)-bin_size+1,bin_size))
    ])
    lst3 = np.array([
        (grp.iloc[:,1].mean(), phred(recall_score(grp.iloc[:,0], grp.iloc[:,2])))
        for grp in (df.iloc[i:i+bin_size] for i in range(0,len(df)-bin_size+1,bin_size))
    ])
    #plt.plot(*lst.T, '--bo')
    #plt.xlabel("Predicted Probability")
    #plt.ylabel("QUAL (Phred-Scaled Recall)")
    return lst, lst2, lst3

# def prob_phred(results):
#     df = pd.read_csv(results, sep="\t", header=0, index_col=["CHROM", "POS"], usecols=['CHROM', 'POS', 'breakca~truth', 'breakca~prob.1', 'breakca~CLASS:']).sort_values(by='breakca~prob.1')
#     df['phred'] = -10*np.log10(1-df['breakca~prob.1'])
#     prob1 = df[df['phred']>-10*np.log10(0.5)]
#     plt.plot(range(len(prob1['phred'])), prob1['phred'])


def strip_type(caller):
    """
        strip the -indel or -snp from the end of a caller name
    """
    if caller.endswith('-snp'):
        return caller[:-len('-snp')]
    elif caller.endswith('-indel'):
        return caller[:-len('-indel')]
    else:
        return caller

def isnan(val):
    return type(val) is float and np.isnan(val)

def get_calls(prepared, callers=None):
    """
        get the alleles in each row of prepared at the location (CHROM/POS) of loc
        when choosing an alt allele, choose from the callers in the order given
    """
    # keep track of the contigs that we've seen
    contigs = set()
    # retrieve the first CHROM/POS location
    loc, predict = yield
    # iterate through each row in the df and check whether they match loc
    for chunk in prepared:
        # if callers is None, retrieve the callers from the columns of the dataframe
        if callers is None:
            callers = [
                caller[:-len('~ALT')] for caller in chunk.columns
                if caller.endswith('~ALT') and not caller.startswith('pg-')
            ]
        if callers is not dict:
            callers = {caller:strip_type(caller) for caller in callers}
        # now iterate through each row (and also convert the POS column to an int)
        for idx, row in chunk.iterrows():
            # check if we already passed the row -- ie we either:
            # 1) we moved onto a new contig or
            # 2) we moved passed the position
            while (
                idx[0] != loc[0] and loc[0] in contigs
            ) or (
                idx[0] == loc[0] and idx[1] > loc[1]
            ):
                # return None if we couldn't find loc in the df
                loc, predict = yield None
            if idx == loc:
                # we found it!
                # now, we must figure out which caller to get the alleles from
                for caller in callers:
                    ref, alt = row[caller+"~REF"], row[caller+"~ALT"]
                    # TODO: make this work for an arbitrary number of variant types for multilabel classification using the other CLASS values in classified
                    # right now, this only works if there's a single binary label
                    if not isnan(ref) and (
                        (isnan(alt) + predict) % 2
                    ):
                        break
                loc, predict = yield callers[caller], ref, alt
            # save the current contig so that we know which ones we've seen
            contigs.add(idx[0])

def main(prepared, classified, callers=None, cs=1000, all_sites=False):
    """
        use the results of the prepare pipeline and the classify pipeline
        to create a VCF with all of the classified sites
    """
    # first, get a generator that can read through each call in the prepared df
    prepared = get_calls(
        pd.read_csv(
            prepared, sep="\t", header=0, index_col=["CHROM", "POS"],
            dtype=str, chunksize=cs, na_values="."
        ), callers
    )
    # flush out the first item in the generator
    next(prepared)
    # also retrieve the classifications as a df
    classified = pd.read_csv(
        classified, sep="\t", header=0, index_col=["CHROM", "POS"],
        dtype={'CHROM':str, 'POS':int}, chunksize=cs, na_values=".",
        usecols=['CHROM', 'POS', 'breakca~truth', 'breakca~prob.1', 'breakca~CLASS:']
    )
    # keep track of how many sites in the classifications df we've had to skip
    skipped = 0
    # iterate through each site in the classifications df, get its alleles, and
    # then return them in a nice-looking dictionary
    for chunk in classified:
        for idx, row in chunk.iterrows():
            try:
                # get the alleles for this CHROM/POS location
                call = prepared.send((idx, row['breakca~CLASS:']))
            except StopIteration:
                call = None
            # check: does the site appear in the prepared pipeline?
            # and does this site have a variant?
            if call is None or (not all_sites and isnan(call[2])):
                skipped += 1
                continue
            caller, ref, alt = call
            # TODO: scale qual by trained linear regression model
            qual = -10*np.log10(1-row['breakca~prob.1'])
            # construct a dictionary with all of the relevant details
            yield {
                'contig': str(idx[0]), 'start': idx[1], 'stop': idx[1]+len(ref),
                'qual': qual, 'alleles': (ref, "." if isnan(alt) else alt), 'info': {'CALLER':caller}
            }
    if skipped:
        warnings.warn(
            "Ignored {:n} classification sites that didn't have a variant.".format(skipped)
        )

def write_vcf(out, records):
    """
        write the records to the output vcf
    """
    vcf = pysam.VariantFile(args.out, mode='w')
    # write the necessary VCF header info
    vcf.header.info.add("CALLER", 1, 'String', "The caller from which this site was taken", Source='varca')
    contigs = set()
    for rec in records:
        try:
            # parse the record into a pysam.VariantRecord
            record = vcf.new_record(
                **rec, samples=None, id=None, filter=None
            )
        except ValueError:
            # add the contig if it hasn't already been added
            if rec['contig'] not in contigs:
                vcf.header.contigs.add(rec['contig'])
                contigs.add(rec['contig'])
            else:
                raise
            # now, try again
            record = vcf.new_record(
                **rec, samples=None, id=None, filter=None
            )
        # write the record to a file
        vcf.write(record)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--out", default=sys.stdout, help="the filename to save the vcf (or bcf) to"
    )
    parser.add_argument(
        "classified", type=Path, help="a sorted, results.tsv.gz file from the output of the classify pipeline"
    )
    parser.add_argument(
        "prepared", type=Path, nargs="?", default=sys.stdin, help="a sorted, merge.tsv.gz file from the prepare pipeline"
    )
    parser.add_argument(
        '-c', "--callers", default=None, help="a comma separated list of the callers from which to choose alleles, supplied in order of priority (default: all of the callers in the file, in the order they appear)"
    )
    parser.add_argument(
        '-s', "--chunksize", type=np.uint32, default=1000, help="how many rows to read into memory at once"
    )
    parser.add_argument(
        '-a', '--all', action='store_true', help="whether to also write non-variant sites to create a gVCF (default: no)"
    )
    args = parser.parse_args()

    write_vcf(args.out, main(args.prepared, args.classified, args.callers, args.chunksize, args.all))
