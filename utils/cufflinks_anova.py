#!/usr/bin/env python
# encoding: utf-8
"""
Perform ANOVA on cufflinks output between two groups, a and b.

    $ python cufflinks_anova.py -v -a sample1 -a sample2 -b sample3 -b sample4 > results.tab

"""
import os
import sys
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
from scipy import stats as st
from toolshed import reader


__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def builddf(a_fpkm, b_fpkm, verbose=False):
    """converts genes.fpkm_tracking files to dataframe."""
    if verbose: sys.stderr.write(">> Building dataframe.\n")
    anames = []
    bnames = []
    genenames = {}
    cases = {}
    for i, f in enumerate(a_fpkm + b_fpkm):
        # assumes genes_fpkm.tracking file is located in a folder named after the sample
        name = os.path.abspath(os.path.join(f, os.path.pardir)).rsplit("/",1)[1]
        if i < len(a_fpkm):
            anames.append(name)
        else:
            bnames.append(name)
        if verbose: sys.stderr.write(">> Adding %s.\n" % name)
        cases[name] = {}
        for c in reader(f):
            if i == 0:
                genenames[c['gene_id']] = c['gene_short_name']
            if c['FPKM_status'] != "OK": continue
            if c['FPKM_conf_lo'] == '0':
                cases[name][c['gene_id']] = 0.
            else:
                cases[name][c['gene_id']] = float(c['FPKM'])
    return anames, bnames, genenames, pd.DataFrame(cases)


def foldchange(a, b):
    amean = np.mean(a)
    bmean = np.mean(b)
    fc = 0.
    if amean != 0. and bmean != 0.:
        if amean > bmean:
            fc = np.power(2, np.log2(amean) - np.log2(bmean)) * -1
        if amean < bmean:
            fc = np.power(2, np.log2(bmean) - np.log2(amean))
    else:
        if amean > bmean:
            fc = amean
        if amean < bmean:
            fc = bmean
    return amean, bmean, fc


def main(args):
    anames, bnames, genenames, df = builddf(args.groupa, args.groupb, args.verbose)
    
    # header with gene and sample info
    header = ['#gene_id', 'common_name']
    header.extend(["a_%s" % name for name in anames])
    header.extend(["b_%s" % name for name in bnames])
    header.extend(['a_mean', 'b_mean', 'foldchange', 'p-value'])
    print args.sep.join(map(str, header))
    
    if args.verbose: sys.stderr.write(">> Performing calculations.\n")

    if args.verbose: 
        sys.stderr.write(">> Correlation table:\n")
        sys.stderr.write("%s\n" % df.corr())
    
    for geneid, values in df.iterrows():
        a = np.array([values[n] for n in anames])
        b = np.array([values[n] for n in bnames])
        fval, pval = st.f_oneway(a, b)
        if str(pval) == 'nan':
            pval = 1.0
        amean, bmean, fc = foldchange(a, b)

        fields = [geneid, genenames.get(geneid)]
        fields.extend([values[n] for n in anames])
        fields.extend([values[n] for n in bnames])
        fields.extend([amean, bmean, fc, pval])
        print args.sep.join(map(str, fields))


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(__doc__, 
            # formatter_class=argparse.ArgumentDefaultsHelpFormatter,
            epilog='Parent directory of genes.fpkm_tracking must be unique. \
            Dependencies are numpy, scipy, pandas, and toolshed.')
    p.add_argument('-a', '--groupa', action='append', required=True,
            help='genes.fpkm_tracking file; use option per file')
    p.add_argument('-b', '--groupb', action='append', required=True,
            help='genes.fpkm_tracking file; use option per file')
    p.add_argument('-s', '--sep', default="\t", 
            help='optional output field separator')
    p.add_argument('-v', '--verbose', default=False, action='store_true', 
            help='for maximum verbosity')
    args = p.parse_args()
    main(args)