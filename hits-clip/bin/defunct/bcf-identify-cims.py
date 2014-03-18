#!/usr/bin/env python
# encoding: utf-8
"""
bcf-identify-cims.py

Created by Joe Brown on 2012-03-19.
"""
import argparse
import os
import sys
from toolshed import nopen
import itertools
import re


def reader(fname, sep="\t"):
    """Reads file, yields dict based on header."""
    line_gen = (l.rstrip("\r\n").split() for l in nopen(fname))
    for line in line_gen:
        header = line
        if header[0].startswith("##"): 
            continue
        break
    header[0] = header[0].lstrip("#")
    for toks in line_gen:
        yield dict(itertools.izip(header, toks))


def identify_cims(vcf):
    """    
    chr start stop read_depth
    """
    variants = []
    for t in reader(vcf):
        start = int(t['POS'])
        stop = start + len(t['ALT'])
        depth = int(re.findall("DP=(\d+)", t['INFO'])[0])
        overlap = False
        #check if current variant overlaps previous, make them not overlap
        if len(variants) > 0 and start >= variants[-1][1] and start <= variants[-1][2]:
            overlap = True
            variants[-1][2] = stop
            previousdepth = variants[-1][3]
            if depth > previousdepth:
                variants[-1][3] = previousdepth + (depth - previousdepth)
        if overlap: continue
        variant = [t['CHROM'], start, stop, depth]
        variants.append(variant)
    for var in variants:
        print "\t".join(str(v) for v in var)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group('required arguments')
    p.add_argument("-i", dest="input", nargs="+", help="Variants call format files.")
    args = p.parse_args()
    if args.input:
        for f in args.input:
            identify_cims(f)
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()