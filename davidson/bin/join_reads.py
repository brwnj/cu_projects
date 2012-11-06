#!/usr/bin/env python
# encoding: utf-8
"""
Join R1 and R2 into SSAKE compatible fasta.
"""
import sys
from cogent.parse.fastq import MinimalFastqParser
from toolshed import nopen

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def fastqtodict(fastq):
    fdict = {}
    for name, seq, qual in MinimalFastqParser(nopen(fastq), strict=False):
        fdict[name] = seq
    return fdict


def main(args):
    r2 = fastqtodict(args.R2)
    for name, seq, qual in MinimalFastqParser(nopen(args.R1), strict=False):
        try:
            r2seq = r2.get(name)
            print ">%s\n%s:%s" % (name, seq, r2seq)
        except KeyError:
            sys.stderr.write(">> No match found for: %s\n" % read)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('R1')
    p.add_argument('R2')
    args = p.parse_args()
    
    main(args)