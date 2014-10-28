#!/usr/bin/env python
# encoding: utf-8
"""
Annotating miRNA bed with features provided by reference bed files and writing
a Cytoscape attributes file.
"""

import argparse
import sys
from pybedtools import BedTool
from itertools import izip

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def main():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('bed', help='bed with miRNA as name')
    p.add_argument('--reference-beds', dest='reference', nargs='+', 
        help='reference beds for each feature to annotate')
    p.add_argument('--names', nargs='+', 
        help='names corresponding to reference files')
    args = p.parse_args()
    if not args.names and not args.reference:
        sys.exit(p.print_help())
    
    bed = BedTool(args.bed)
    
    # create the reference beds
    reference = {}
    for refname, refbed in izip(args.names, args.reference):
        reference[refname] = BedTool(refbed)
    
    for refname in args.names:
    
        # intersect the mirna bed with the reference annotations
        for b in bed.intersect(reference[refname], s=True, stream=True):
            # Cytoscape formatting
            fields = (b.name, "=", refname)
            print " ".join(map(str, fields))


if __name__ == "__main__":
    main()