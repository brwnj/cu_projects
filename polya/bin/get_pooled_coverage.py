#!/usr/bin/env python
# encoding: utf-8
"""
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
from bsub import bsub
from toolshed import reader
from itertools import combinations
from collections import defaultdict
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

def get_sample_name(fname, pattern):
    """
    >>>get_sample_name("PK96.pos.bedgraph.gz", ".bedgraph")
    'PK96.pos'
    """
    return fname.split(pattern)[0]

def get_compression_setting(fname):
    return "gzip" if fname.endswith(".gz") else None

def _sf_deseq(counts):
    """
    Calculate DESeq scaling factor per sample.
    """
    # masked array to discard inf, -inf, and nan
    ma = np.ma.masked_invalid(counts)
    return np.exp(np.ma.median(ma))

def norm_deseq(df):
    """
    Normalize by DESeq scaling factor, which is computed as the median of
    the ratio, for each row (gene), of its read count over its geometric
    mean across all samples. Return new counts dataframe.

    Details:
    --------
    http://genomebiology.com/2010/11/10/R106

    Parameters:
    -----------
    - df: pandas dataframe.
    """
    # log of counts
    lg = df.apply(np.log)
    # per sample: exponential(median(log(counts) - geometric mean))
    sf = lg.sub(lg.mean(axis=1), axis=0).apply(_sf_deseq, axis=0)
    # apply scaling
    df = df.div(sf, axis=1)
    return df

def main(bedgraphs, metadata):
    pools = defaultdict(list)
    for toks in reader(metadata):
        for k, v in toks.iteritems():
            if k.startswith("Pool") and v == "TRUE":
                # get the samples
                pool_name = k.split("_")[-1]
                pools[pool_name].append(toks['alias'])
    for pool, samples in pools.iteritems():
        for strand in ["pos", "neg"]:
            print >>sys.stderr, ">> processing", pool
            files = [f for f in bedgraphs if os.path.basename(f).split(".")[0] in samples and strand in os.path.basename(f) and "UMI" not in f]
            if len(files) == 0: continue
            assert len(files) == len(samples), "All count files not found for {pool}".format(pool=pool)
            df_list = [pd.read_table(f, names=["chrom", "start", "stop",get_sample_name(f, ".bedgraph")], index_col=["chrom", "start", "stop"], compression="gzip") for f in files]
            combined_df = None
            for i, (dfa, dfb) in enumerate(combinations(df_list, 2)):
                if i > len(df_list) / 2: continue
                if i == 0:
                    combined_df = dfa
                    combined_df = combined_df.join(dfb)
                else:
                    combined_df = combined_df.join(dfb)
            # normalize the counts
            combined_df = norm_deseq(combined_df)
            # round the normalized counts up to int
            # don't want to throw out single counts at any site
            combined_df = combined_df.apply(np.ceil)
            combined_df.fillna(0, inplace=True)
            # sum the rows
            combined_df[pool] = combined_df.sum(axis=1)
            # print results
            combined_df[pool].astype('int').to_csv("{pool}.{strand}.bedgraph".format(pool=pool, strand=strand), sep="\t")

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("metadata")
    p.add_argument("bedgraphs", nargs="+", help="sorted sample bedgraph files")
    args = p.parse_args()
    main(args.bedgraphs, args.metadata)
