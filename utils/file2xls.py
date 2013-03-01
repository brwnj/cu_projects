#!/usr/bin/env python
# encoding: utf-8

import os.path as op
import pandas as pd
from toolshed import nopen

def main(args):
    for f in args.txt:
        fin = nopen(f)
        df = pd.io.parsers.read_table(fin, header=args.header, verbose=args.verbose)
        df.to_excel(op.splitext(f)[0] + ".xls", index=args.index, 
                    na_rep=args.narep)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('txt', nargs="+", help='text file')
    p.add_argument('--header', type=int, default=None,
            help='specify header row [ %(default)s ]')
    p.add_argument('--index', action='store_true',
            help='print row names [ %(default)s ]')
    p.add_argument('--narep', default="",
            help='missing data representation [ %(default)s ]')
    p.add_argument('--verbose', action='store_true',
            help='maximum verbosity [ %(default)s ]')
    main(p.parse_args())