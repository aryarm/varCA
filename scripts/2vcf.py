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



def plot_line(lst):
    roc = lst.T
    p = np.poly1d(np.polyfit(roc[0], roc[1], 1))
    plt.scatter(roc[0], roc[1], label="_nolegend_")
    plt.xlabel("Probability")
    plt.ylabel("TPR")
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

def plot_tpr_probs(results, bin_size=20):
    """ bin the sites and calculate an accuracy for that bin """
    df = pd.read_csv(results, sep="\t", header=0, index_col=["CHROM", "POS"], usecols=['CHROM', 'POS', 'breakca~truth', 'breakca~prob.1', 'breakca~CLASS:']).sort_values(by='breakca~prob.1')
    # split the dataset into high (above 0.5) vs low (below 0.5) scores
    df1 = df[df['breakca~truth']>=0.5]
    df2 = df[df['breakca~truth']<0.5][::-1]
    lst = np.array([
        (grp.iloc[:,1].mean(), recall_score(grp.iloc[:,0], grp.iloc[:,2]))
        for grp in (df1.iloc[i:i+bin_size] for i in range(0,len(df1)-bin_size+1,bin_size))
    ])
    bin_size = int(bin_size*len(df2)/len(df1))
    # lst2 = np.array([
    #     (grp.iloc[:,1].mean(), phred(recall_score(grp.iloc[:,0], grp.iloc[:,2], pos_label=0)))
    #     for grp in (df2.iloc[i:i+bin_size] for i in range(0,len(df2)-bin_size+1,bin_size))
    # ])
    # lst3 = np.array([
    #     (grp.iloc[:,1].mean(), phred(recall_score(grp.iloc[:,0], grp.iloc[:,2])))
    #     for grp in (df.iloc[i:i+bin_size] for i in range(0,len(df)-bin_size+1,bin_size))
    # ])
    #plt.plot(*lst.T, '--bo')
    #plt.xlabel("Predicted Probability")
    #plt.ylabel("QUAL (Phred-Scaled Recall)")
    return lst
    return lst, lst2, lst3

# def plot_tpr_probs(results, bin_size):
#     """ bin the sites and calculate an accuracy for that bin """
#     df = pd.read_csv(results, sep="\t", header=0, index_col=["CHROM", "POS"], usecols=['CHROM', 'POS', 'breakca~truth', 'breakca~prob.1', 'breakca~CLASS:']).sort_values(by='breakca~prob.1')
#     # split the dataset into high (above 0.5) vs low (below 0.5) scores
#     df1 = df[df['breakca~truth']>=0.5]
#     df2 = df[df['breakca~truth']<0.5][::-1]
#     lst = np.array([
#         (grp[1].iloc[:,1].mean(), recall_score(grp[1].iloc[:,0], grp[1].iloc[:,2]), len(grp[1]))
#         for grp in df1.groupby(pd.cut(df1['breakca~prob.1'], pd.interval_range(0.5, 1, bin_size)))
#     ])
#     bin_size = int(bin_size*len(df2)/len(df1))
#     # lst2 = np.array([
#     #     (grp[1].iloc[:,1].mean(), recall_score(grp[1].iloc[:,0], grp[1].iloc[:,2], pos_label=0))
#     #     for grp in df1.groupby(pd.cut(df1['breakca~prob.1'], pd.interval_range(0, 0.5, bin_size)))
#     # ])
#     #plt.plot(*lst.T, '--bo')
#     #plt.xlabel("Predicted Probability")
#     #plt.ylabel("QUAL (Phred-Scaled Recall)")
#     return lst

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
        caller = caller[:-len('-snp')]
    elif caller.endswith('-indel'):
        caller = caller[:-len('-indel')]
    # if there is still a dash, get everything after it
    i = caller.rfind('-')
    if i != -1:
        caller = caller[i+1:]
    return caller

def isnan(val):
    return type(val) is float and np.isnan(val)

def get_calls(prepared, callers=None, pretty=False):
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
            callers = {caller: (strip_type(caller) if pretty else caller) for caller in callers}
        # now iterate through each row (and also convert the POS column to an int)
        for idx, row in chunk.iterrows():
            # check if we already passed the row -- ie we either:
            # 1) moved onto a new contig or
            # 2) moved passed the position
            while (
                idx[0] != loc[0] and loc[0] in contigs
            ) or (
                idx[0] == loc[0] and idx[1] > loc[1]
            ):
                # return None if we couldn't find loc in the df
                loc, predict = yield None
            if idx == loc:
                # we found it!
                found = False
                # now, we must figure out which caller to get the alleles from
                for caller in callers:
                    ref, alt = row[caller+"~REF"], row[caller+"~ALT"]
                    # TODO: make this work for an arbitrary number of variant types for multilabel classification using the other CLASS values in classified
                    # right now, this only works if there's a single binary label
                    if not isnan(ref) and (
                        (isnan(alt) + predict) % 2
                    ):
                        found = True
                        break
                if found:
                    loc, predict = yield callers[caller], ref, alt
                else:
                    # if we know there is a variant here, but none of the other
                    # callers found it, just label it as a non-variant
                    # TODO: figure out the alt allele from inspecting the ref genome?
                    loc, predict = yield 'varca', row['REF'], float('nan')
            # save the current contig so that we know which ones we've seen
            contigs.add(idx[0])

def prob2qual(prob):
    # TODO: scale qual by trained linear regression model
    return -10*np.log10(1-prob)

def main(prepared, classified, callers=None, cs=1000, all_sites=False, pretty=False):
    """
        use the results of the prepare pipeline and the classify pipeline
        to create a VCF with all of the classified sites
    """
    # first, get a generator that can read through each call in the prepared df
    prepared = get_calls(
        pd.read_csv(
            prepared, sep="\t", header=0, index_col=["CHROM", "POS"],
            dtype=str, chunksize=cs, na_values="."
        ), callers, pretty
    )
    # flush out the first item in the generator
    next(prepared)
    # also retrieve the classifications as a df
    classified = pd.read_csv(
        classified, sep="\t", header=0, index_col=["CHROM", "POS"],
        dtype={'CHROM':str, 'POS':int}, chunksize=cs, na_values="."
    )
    # keep track of how many sites in the classifications df we've had to skip
    skipped = 0
    # keep track of how many sites we skipped but were predicted to have a variant
    no_alts = 0
    # iterate through each site in the classifications df, get its alleles, and
    # then return them in a nice-looking dictionary
    for chunk in classified:
        for idx, row in chunk.iterrows():
            try:
                # get the alleles for this CHROM/POS location
                call = prepared.send((idx, row['breakca~CLASS:']))
            except StopIteration:
                call = None
            # we found a variant but couldn't find an alternate allele!
            no_alts += not (call is None or (isnan(call[2]) + row['breakca~CLASS:']) % 2)
            # check: does the site appear in the prepared pipeline?
            # and does this site have a variant?
            if call is None or (not all_sites and isnan(call[2])):
                skipped += 1
                continue
            caller, ref, alt = call
            qual = prob2qual(row['breakca~prob.'+str(int(not isnan(alt)))])
            # construct a dictionary with all of the relevant details
            yield {
                'contig': str(idx[0]), 'start': idx[1], 'stop': idx[1]+len(ref),
                'qual': qual, 'alleles': (ref, "." if isnan(alt) else alt), 'info': {'CALLER':caller}
            }
    if skipped:
        warnings.warn(
            "Ignored {:n} classification sites that didn't have a variant.".format(skipped)
        )
    if no_alts:
        warnings.warn(
            "Ignored {:n} sites that we predicted to have variants but didn't appear in any of the callers.".format(no_alts)
        )

def write_vcf(out, records):
    """
        write the records to the output vcf
    """
    vcf = pysam.VariantFile(args.out, mode='w')
    # write the necessary VCF header info
    vcf.header.info.add("CALLER", 1, 'String', "The caller from which this site was taken")
    contigs = set()
    for rec in records:
        # handle pysam increasing the start and end sites by 1
        rec['start'] -= 1
        rec['stop'] -= 1
        # parse the record into a pysam.VariantRecord
        try:
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
    return vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--out", default=sys.stdout, help="the filename to save the vcf (or bcf) to"
    )
    parser.add_argument(
        "classified", type=Path, help="a sorted, results.tsv.gz file from the output of the classify pipeline"
    )
    parser.add_argument(
        "prepared", type=Path, nargs="?", default=sys.stdin, help="a sorted, merge.tsv.gz file from the prepare pipeline (if not supplied, this is read from stdin)"
    )
    parser.add_argument(
        '-c', "--callers", default=None, help="a comma separated list of the callers from which to choose alleles, supplied in order of priority (default: all of the callers in the file, in the order they appear)"
    )
    parser.add_argument(
        '-s', "--chunksize", type=np.uint32, default=100000, help="how many rows to read into memory at once (default: 100,000)"
    )
    parser.add_argument(
        '-a', '--all', action='store_true', help="whether to also write non-variant sites to create a gVCF (default: no)"
    )
    parser.add_argument(
        '-p', '--pretty', action='store_true', help="should caller names appear in the vcf by their pretty form (with all dashes intelligently removed) or their original caller ID form? (default: the original form)"
    )
    parser.add_argument(
        '-n', '--nothing', action='store_true', help="do absolutely nothing, except read the arguments (for testing and internal use)"
    )
    args = parser.parse_args()

    callers = None
    if args.callers is not None:
        callers = args.callers.split(",")

    if not args.nothing:
        vcf = write_vcf(args.out, main(args.prepared, args.classified, callers, args.chunksize, args.all, args.pretty))
    else:
        plt.ion()
        lst = plot_tpr_probs(args.classified)
