#!/usr/bin/env python
# encoding: utf-8
"""
takes gene_exp.diff files and generates matrix of fold change by significant
gene.

updating this script to transform genes.fpkm_tracking to data matrix
"""
import sys
from toolshed import reader
from collections import defaultdict

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def main(args):
    gm = defaultdict(dict)
    for trackingfile in args.trackingfiles:
        # reading expression profiles
        for gene in reader(trackingfile):
            #if e['significant'] == "yes" and float(e['p_value']) < 0.05:
            gm[trackingfile.split("_")[0]][gene['gene_id']] = {}
            gm[trackingfile.split("_")[0]][gene['gene_id']] = gene[args.column]

    # print the matrix
    caselist = sorted(gm.keys())
    genes = sorted(gm[caselist[0]].keys())
    print "#gene_id\t" + "\t".join(k for k in caselist)
    for gene in genes:
        print gene + "\t" + "\t".join(gm[c][gene] for c in caselist)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('trackingfiles', nargs='+', 
            help='space delimited list or glob of genes.fpkm_tracking files')
    # p.add_argument('--pvalue', '-p', default=0.05, help='p-value cutoff')
    p.add_argument('--column', '-c', default='FPKM', 
            help='column to join')
    args = p.parse_args()
    main(args)