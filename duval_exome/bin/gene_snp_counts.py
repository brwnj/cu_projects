#!/usr/bin/env python
# coding=utf-8
"""
"""

import pandas as pd
import os.path as op
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def xreftodict(xref):
    d = {}
    for t in reader(xref, header=["ens","sym"]):
        d[t['ens']] = t['sym']
    return d


def main(snps, xref, coverage):
    xr = xreftodict(xref)
    filename, ext = op.splitext(snps)
    df = pd.read_table(snps, compression="gzip" if ext == ".gz" else None, skiprows=2)
    df = df[df.Effect == "NON_SYNONYMOUS_CODING"]
    df = df[df.Coverage >= coverage]
    for gene, gdf in df.groupby('Gene_ID'):
        try:
            sym = xr[gene]
        except KeyError:
            sym = gene
        print "\t".join([gene, sym, str(len(gdf))])


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('snps', help="snpeff annotated txt file")
    p.add_argument('xref', help="ens gene id [tab] gene symbol")
    p.add_argument('--coverage', default=10, type=int, help="minimum coverage to be considered")
    args = vars(p.parse_args())
    main(**args)
