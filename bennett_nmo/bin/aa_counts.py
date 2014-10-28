#!/usr/bin/env python
# encoding: utf-8
"""
Amino Acid sequence counts of unique CDR3 sequences parsed from IMGT
AA-sequences.
"""
import os
import sys
import gzip
import pandas as pd
from collections import Counter
from toolshed import reader, nopen
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def gsample(fname):
    # incoming name format:
    # 5_AA-sequences_ON10_03B_gc1_170114.txt
    # 5_AA-sequences_control_ACATCG_170114.txt
    f = os.path.basename(fname)
    return f.split("_", 2)[-1].rsplit("_", 1)[0]


def main(aaseqs):
    # count the unique sequences
    counters = {}
    for f in aaseqs:
        sample = gsample(f)
        print >>sys.stderr, ">> reading in", f
        counters[sample] = Counter()
        for l in reader(f):
            if l['Functionality'] != "productive": continue
            rname, umi, cregion, fwork = l['Sequence ID'].split(':')
            cdr3 = l['CDR3-IMGT']
            counters[sample].update(["{cdr3}:{cregion}_{fwork}".format(**locals())])
    # prep to create dataframe
    predf = {}
    for name, counts in counters.iteritems():
        predf[name] = {}
        for seq, count in counts.iteritems():
            predf[name][seq] = count
    # create dataframe
    df = pd.DataFrame(predf)
    # create multiindex via split
    df.index = pd.MultiIndex.from_tuples([x.split(":") for x in df.index], names=['CDR3','C_Fwork'])
    df.to_csv(sys.stdout, sep="\t", na_rep=0)


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__,
            formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument("aaseqs", nargs="+", help="IMGT AA-sequences")
    args = p.parse_args()
    main(args.aaseqs)
