#!/usr/bin/env python
# encoding: utf-8
"""
Amino Acid sequence counts of unique CDR3 sequences parsed from IMGT
AA-sequences.
"""
import os
import csv
import sys
import gzip
import argparse
import pandas as pd
from itertools import izip
from collections import Counter

def nopen(f, mode="rb"):
    f = os.path.expanduser(os.path.expandvars(f))
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
         else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
         else open(f, mode)

def reader(fname, header=True, sep="\t"):
    dialect = csv.excel
    dialect.delimiter = sep
    line_gen = csv.reader(nopen(fname), dialect=dialect)
    a_dict = dict
    header = line_gen.next()
    header[0] = header[0].lstrip("#")
    for toks in line_gen:
        yield a_dict(izip(header, toks))

def gsample(fname):
    # incoming name format: 5_AA-sequences_2_ACATCG_1_270613.txt
    fname = os.path.basename(fname)
    lst = fname.split("_")
    patient = lst[2]
    barcode = lst[3]
    return "{patient}_{barcode}".format(**locals())

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
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("aaseqs", nargs="+", help="IMGT AA-sequences")
    args = p.parse_args()
    main(args.aaseqs)
