#!/usr/bin/env python
# encoding: utf-8
"""
Convert scaled miRNA abundance files into a table by case.
"""
import argparse
import sys
from collections import defaultdict
from toolshed import reader
from pybedtools import BedTool

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def get_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('files', metavar='ABUNDANCE_FILES', nargs="+", help='mirna name and scaled abundance')
    p.add_argument('--reference', '-r', help='mirbase bed with all miRNAs')
    args = p.parse_args()
    return args


def main():
    args = get_args()
    cases = defaultdict(dict)
    for f in args.files:
        # not all files contain all miRNA, so fill in with 0s first
        for b in BedTool(args.reference):
            cases[f[:-4]][b.name] = {}
            cases[f[:-4]][b.name] = '0'
        for mirna in reader(f, header="name abundance".split()):
            cases[f[:-4]][mirna['name']] = mirna['abundance']
            
    caselist = sorted(cases.keys())
    mirnas = [b.name for b in BedTool(args.reference)]
    print "#mirna_id\t" + "\t".join(k for k in caselist)
    for mirna in mirnas:
        print mirna + "\t" + "\t".join(cases[k][mirna] for k in caselist)


if __name__ == "__main__":
    main()