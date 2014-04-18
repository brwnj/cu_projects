#!/usr/bin/env python
# coding=utf-8
"""
"""
import sys
import pandas as pd
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def main(snpeff, effect, min_coverage=10):
    df = pd.read_table(snpeff, header=2, compression="gzip")
    df = df[(df['Effect'].isin(effect)) & (df['Coverage'] >= min_coverage)]
    try:
        df.to_csv(sys.stdout, index=False, na_rep="NA")
    except IOError:
        # allow piping to head without issues
        pass


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('snpeff')
    p.add_argument('-e', '--effect', action="append", default=['NON_SYNONYMOUS_CODING'], help="mutation effect defined by snpeff; can be specified more than once")
    p.add_argument('--min-coverage', default=10, type=int, help="minimum allowable coverage")
    args = vars(p.parse_args())
    main(**args)
