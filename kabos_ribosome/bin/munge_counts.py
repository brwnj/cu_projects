#!/usr/bin/env python
# coding=utf-8
"""
"""
import sys
import pandas as pd
import os.path as op
from argparse import ArgumentParser


def main(counts):
    df = pd.read_table(counts, header=1)
    # drop unnecessary columns
    df = df.drop("Chr", 1)
    df = df.drop("Start", 1)
    df = df.drop("End", 1)
    df = df.drop("Strand", 1)
    df = df.drop("Length", 1)
    # remove name of first column
    df.rename(columns={'Geneid': ''}, inplace=True)
    # rename sample columns to basename minus '.bam'
    df.rename(columns=lambda x: op.basename(x).replace("-", "_").rstrip(".bam"), inplace=True)
    try:
        df.to_csv(sys.stdout, sep="\t", na_rep=0, index=False)
    except IOError:
        pass


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__)
    p.add_argument('counts', help="featureCounts count table")
    args = vars(p.parse_args())
    main(**args)
