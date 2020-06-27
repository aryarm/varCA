#!/usr/bin/env python3
import sys
import pysam
import argparse
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
from sklearn.metrics import accuracy_score, precision_score
from sklearn.metrics import roc_curve



def phred(vals):
    """ apply the phred scale to the vals provided """
    return -10*np.log10(1-vals)
    return -10*np.ma.log10(1-vals).filled(-3) # fill all infinite values with a phred scale of 30

def plot_line(lst, show_discards=False):
    plt.clf()
    roc = np.copy(lst.T)
    roc[1] = phred(roc[1])
    max_val = phred(max(roc[2])/(max(roc[2])+1))
    # discard inf or na cols
    inf_or_na_cols = np.isinf(roc).any(axis=0) | np.isnan(roc).any(axis=0)
    # warn the user if we're discarding the majority of points
    discarded = np.sum(inf_or_na_cols)/roc.shape[1]*100
    if (not show_discards and discarded != 0) or discarded >= 50:
        warnings.warn("Discarding NaN or Inf points ({}% of points)".format(discarded))
    roc = roc[:,~(inf_or_na_cols)]
    # perform a simple linear regression
    p = np.polyfit(roc[0], roc[1], 1)
    r_squared = 1 - (sum((roc[1] - (p[0] * roc[0] + p[1]))**2) / ((len(roc[1]) - 1) * np.var(roc[1], ddof=1)))
    p = np.poly1d(p)
    # plot the points and the line
    plt.scatter(roc[0], roc[1], color='r', label="_nolegend_")
    if max(roc[0]) <= 1:
        plt.xlabel("RF Probability")
    elif max(roc[0]) <= np.pi/2:
        plt.xlabel("Reverse Arcsin of RF Probability")
    else:
        plt.xlabel("Phred-Scaled RF Probability")
    plt.ylabel("Phred-Scaled Accuracy (QUAL)")
    plt.plot(
        roc[0],
        p(roc[0]),
        label=str(p)+"\nr-squared: "+str(round(r_squared, 2))+ \
            ("\ndiscarded: "+str(int(discarded))+"%" if show_discards else "")
    )
    plt.hlines(max_val, min(roc[0]), max(roc[0]), colors='g', linestyles='dashed')
    plt.legend(frameon=False, loc='lower right')
    plt.tight_layout()

def eqbin_mean(grp, log=True, pseudo=True, discards_ok=False, inverse=False):
    if inverse:
        return np.arcsin(grp.mean())
    else:
        if log:
            if discards_ok or not pseudo:
                return phred(grp).mean()
            else:
                return phred(grp.sum()/(len(grp) + pseudo))
        else:
           return grp.mean()

def tpr_probs(df, bins=15, eqbin=True, log=True, pseudo=True, discards_ok=False, inverse=False):
    """ bin the sites and calculate an accuracy (predicted positive value) for that bin """
    # retrieve only predicted positives
    df = df[df['varca~CLASS:']>=0.5]
    if eqbin:
        bin_size = int(len(df)/bins)
        # create bins (ie each grp) and add a single false positive to it so we don't get Inf
        lst = np.array([
            (
                eqbin_mean(grp['varca~prob.1'], log, pseudo, discards_ok, inverse),
                precision_score(
                    np.append(grp['varca~truth'].values, 0) if pseudo else grp['varca~truth'],
                    np.append(grp['varca~CLASS:'].values, 1) if pseudo else grp['varca~CLASS:']
                ),
                grp['varca~prob.1'].size
            )
            for grp in (df.iloc[i:i+bin_size] for i in range(0,len(df)-bin_size+1,bin_size))
        ])
    else:
        if log:
            df = df.copy()
            df['varca~prob.1'] = phred(df['varca~prob.1'])
        start = phred(0.5) if log else 0.5
        # get the end excluding inf values (in case log == True)
        end = df.loc[df['varca~prob.1'] != np.inf, 'varca~prob.1'].max()
        # create bins (ie each grp) and add a single false positive to it so we don't get Inf
        lst = np.array([
            (
                grp[1]['varca~prob.1'].mean(),
                precision_score(
                    np.append(grp[1]['varca~truth'].values, 0) if pseudo else grp['varca~truth'],
                    np.append(grp[1]['varca~CLASS:'].values, 1) if pseudo else grp['varca~CLASS:']
                ),
                grp[1]['varca~prob.1'].size
            )
            for grp in df.groupby(pd.cut(df['varca~prob.1'], pd.interval_range(start, end, bins)))
        ])
    return lst



def strip_type(caller):
    """
        strip the -indel or -snp from the end of a caller name
    """
    vartype = ''
    if caller.endswith('-snp'):
        caller = caller[:-len('-snp')]
        vartype = 'snp'
    elif caller.endswith('-indel'):
        caller = caller[:-len('-indel')]
        vartype = 'indel'
    # if there is still a dash, get everything after it
    i = caller.rfind('-')
    if i != -1:
        caller = caller[i+1:]
    return caller, vartype

def isnan(val):
    return type(val) is float and np.isnan(val)

def get_calls(prepared, callers=None, pretty=False):
    """
        get the alleles in each row of prepared at the location (CHROM/POS) of loc
        when choosing an alt allele, choose from the callers in the order given
    """
    # keep track of the contigs that we've seen
    contigs = set()
    # whether we've read the header yet
    read_header = False
    # iterate through each row in the df and check whether they match loc
    for chunk in prepared:
        # do some special stuff (parsing the header) on the very first iteration
        if not read_header:
            # if callers is None, retrieve the callers from the columns of the dataframe
            if callers is None:
                callers = [
                    caller[:-len('~ALT')] for caller in chunk.columns
                    if caller.endswith('~ALT') and not caller.startswith('pg-')
                ]
            # what types of variants are we dealing with? let's count how many
            # times they appear in the caller names
            vartypes = {'snp': 0, 'indel': 0}
            # also, let's retrieve the callers as a dictionary
            pretty_callers = {}
            for caller in callers:
                pretty_caller, vartype = strip_type(caller)
                # don't beautify the callers if the user didn't request it
                pretty_callers[caller] = pretty_caller if pretty else caller
                # keep track of how often each vartype appears
                if vartype in vartypes:
                    vartypes[vartype] += 1
            callers = pretty_callers
            # retrieve the first CHROM/POS location and yield whether we are reading indels or snps
            loc, predict = yield max(vartypes, key=vartypes.get)
            read_header = True
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

def prob2qual(prob, vartype):
    # values are from linear model that we created from using the "-i" option
    if vartype == 'snp':
        return 0.6237*phred(prob)+8.075
    elif vartype == 'indel':
        return 0.8463*phred(prob)+2.724
    else:
        # we shouldn't ever encounter this situation, but just in case...
        return phred(prob)

def main(prepared, classified, callers=None, cs=1000, all_sites=False, pretty=False, vartype=None):
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
    # flush out the first item in the generator: the vartype
    if vartype is None:
        vartype = next(prepared)
    else:
        # if the user already gave us the vartype, then just discard this
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
                call = prepared.send((idx, row['varca~CLASS:']))
            except StopIteration:
                call = None
            # we found a variant but couldn't find an alternate allele!
            no_alts += not (call is None or (isnan(call[2]) + row['varca~CLASS:']) % 2)
            # check: does the site appear in the prepared pipeline?
            # and does this site have a variant?
            if call is None or (not all_sites and isnan(call[2])):
                skipped += 1
                continue
            caller, ref, alt = call
            qual = prob2qual(
                row['varca~prob.'+str(int(not isnan(alt)))], vartype
            )
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
        '-t', '--type', choices=['indel', 'snp'], default=None, help="whether to recalibrate QUAL values assuming your data are SNPs or indels (default: infer from callers)"
    )
    parser.add_argument(
        '-i', '--internal', action='store_true', help="For testing and internal use: recalibrate the QUAL scores (assumes varca~truth column exists in classified)"
    )
    args = parser.parse_args()

    callers = None
    if args.callers is not None:
        callers = args.callers.split(",")

    if not args.internal:
        vcf = write_vcf(args.out, main(args.prepared, args.classified, callers, args.chunksize, args.all, args.pretty, args.type))
    else:
        if not sys.flags.interactive:
            sys.exit("ERROR: You must run this script in python's interactive mode (using python's -i flag) when providing the -i flag to this script.")
        try:
            df = pd.read_csv(args.classified, sep="\t", header=0, index_col=["CHROM", "POS"], usecols=['CHROM', 'POS', 'varca~truth', 'varca~prob.1', 'varca~CLASS:'], low_memory=False).sort_values(by='varca~prob.1')
        except ValueError:
            df = pd.read_csv(args.classified, sep="\t", header=0, index_col=["CHROM", "POS"], usecols=['CHROM', 'POS', 'breakca~truth', 'breakca~prob.1', 'breakca~CLASS:'], low_memory=False).sort_values(by='breakca~prob.1')
            df.columns = ['varca~truth', 'varca~prob.1', 'varca~CLASS:']
        plt.ion()
        plot_line(tpr_probs(df))
