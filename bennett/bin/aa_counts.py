#!/usr/bin/env python
# encoding: utf-8
"""
Amino Acid sequence counts of unique CDR3 sequences parsed from IMGT
AA-sequences.
"""
import os
import sys
import argparse
import pandas as pd
from toolshed import reader
from collections import Counter, defaultdict

def gsample(fname):
    # incoming name format: 5_AA-sequences_2_ACATCG_1_270613.txt
    fname = os.path.basename(fname)
    lst = fname.split("_")
    patient = lst[2]
    barcode = lst[3]
    return "{patient}_{barcode}".format(**locals())

def main(aaseqs):
    # read in all of the files
    sdata = defaultdict(list)
    for f in aaseqs:
        sample = gsample(f)
        print >>sys.stderr, ">> reading in", f
        for l in reader(f):
            sdata[sample].append(l)
    # count the unique sequences
    counters = {}
    for name, lines in sdata.iteritems():
        counters[name] = Counter()
        for l in lines:
            if l['Functionality'] != "productive": continue
            rname, umi, cregion, fwork = l['Sequence ID'].split(':')
            cdr3 = l['CDR3-IMGT']
            counters[name].update(["{cdr3}:{cregion}_{fwork}".format(**locals())])
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
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("aaseqs", nargs="+", help="IMGT AA-sequences")
    args = p.parse_args()
    main(args.aaseqs)
