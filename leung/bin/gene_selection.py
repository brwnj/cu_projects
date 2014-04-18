#!/usr/bin/env python
# coding=utf-8
"""
"""
import sys
import pandas as pd
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def main(files, effect, min_coverage=10, gene=None, split_sep="_"):
    tables = []

    for f in files:
        sample = f.split(split_sep)[0]
        df = pd.read_table(f, header=2, compression="gzip")
        df['Sample'] = sample
        tables.append(df)

    df = pd.concat(tables)

    # filter
    df = df[(df['Effect'] in effect)\
                & (df['Coverage'] >= min_coverage)]
    df.dropna(how="all", inplace=True)
    if gene:
        # print only selected genes
        df.set_index("Gene_ID", inplace=True)
        df = df.ix[gene]
        df.dropna(how="all", inplace=True)
        df.to_csv(sys.stdout, index=False, na_rep="NA")

    else:
        # print all
        df.to_csv(sys.stdout, index=False, na_rep="NA")


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('files', nargs="+")
    p.add_argument('-e', '--effect', nargs="append", default=['NON_SYNONYMOUS_CODING'], help="mutation effect defined by snpeff; can be specified more than once")
    p.add_argument('-g', '--gene', action="append", help="genes to output; empty for all")
    p.add_argument('--min-coverage', default=10, type=int, help="minimum allowable coverage")
    p.add_argument('--split-sep', help="split pattern on file name to grab sample name")
    # p.add_argument('-')
    args = vars(p.parse_args())
    main(**args)
