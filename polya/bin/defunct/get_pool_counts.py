#!/usr/bin/env python
# encoding: utf-8
"""
Create pooled count files for testing.
"""
import os
import sys
import gzip
import argparse
import numpy as np
import pandas as pd
from itertools import izip
from toolshed import reader
from collections import defaultdict

def get_sample_name(fname):
    # clip off ".counts.txt.gz"
    return fname[:-14]

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

def main(count_files, metadata):
    pools = defaultdict(list)
    for toks in reader(metadata):
        for k, v in toks.iteritems():
            if k.startswith("Pool") and v == "TRUE":
                # get the samples
                pool_name = k.split("_")[-1]
                pools[pool_name].append(toks['alias'])
    for pool, samples in pools.iteritems():
        print >>sys.stderr, ">> processing", pool
        for strand in ["pos", "neg"]:
            files = [f for f in count_files if os.path.basename(f).split(".")[0] in samples and strand in os.path.basename(f)]
            # simplest way to join files into a dataframe
            raw_count_data = {}
            for file_path in files:
                sample = get_sample_name(file_path)
                raw_count_data[sample] = {}
                for toks in reader(file_path, header=['gene', 'site', 'count']):
                    raw_count_data[sample]["{gene}:{site}".format(gene=toks['gene'], site=toks['site'])] = int(toks['count'])
            # dataframe from dict of dicts
            count_data = pd.DataFrame(raw_count_data)
            # will need to split into multiindex here to match new count fmt
            count_data.index = pd.MultiIndex.from_tuples([x.split(":") for x in count_data.index], names=['gene','site'])
            # normalize the counts
            count_data = norm_deseq(count_data)
            # round the normalized counts up to int
            # don't want to throw out single counts at any site
            count_data = count_data.apply(np.ceil)
            # sum the rows
            count_data[pool] = count_data.sum(axis=1)
            # print results
            out_file = gzip.open("{pool}.{strand}.txt.gz".format(pool=pool, strand=strand), "wb")
            count_data[pool].astype('int').to_csv(out_file, sep="\t")
            out_file.close()

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("metadata")
    p.add_argument("counts", nargs="+")
    args = p.parse_args()
    main(args.counts, args.metadata)
