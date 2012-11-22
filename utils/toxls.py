#!/usr/bin/env python
# encoding: utf-8

import os.path as op
import pandas as pd

def main(args):
    df = pd.io.parsers.read_table(args.txt, header=args.header, verbose=args.verbose)
    df.to_excel(op.splitext(args.txt)[0] + ".xls")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('txt', help='text file')
    p.add_argument('--header', type=int, default=None, help='specify header row [ None ]')
    p.add_argument('--verbose', action='store_true', help='maximum verbosity [ False ]')
    # comment indicator
    # delimiter
    # na values
    args = p.parse_args()
    main(args)