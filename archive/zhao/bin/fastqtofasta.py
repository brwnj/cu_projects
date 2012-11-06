#!/usr/bin/env python
# encoding: utf-8
"""
process-polya.py

Created by Joe Brown on 2012-03-13.
"""
import argparse
import os
import sys
from toolshed import nopen
from cogent.parse.fastq import MinimalFastqParser


def fastqtofasta(fastq):
    """Reads fastq and prints fasta to stdout.
    
    """
    for label, seq, qual in MinimalFastqParser(nopen(fastq), strict=False):
        #print fasta
        print '>' + label
        print seq


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", dest="fastq", help="Reads in fastq format. Accepts stdin, .gz, .bz?, urls.")
    args = p.parse_args()
    if args.fastq:
        fastqtofasta(args.fastq)
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()