#!/usr/bin/env python
# encoding: utf-8
"""
Converts gene names to gene symbols with provided cross-reference file.
"""

import argparse
import os
import sys
from toolshed import nopen, reader
from pybedtools import BedTool

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def cross_ref(kgxref, table_id, table_symbol):
    """Returns dictionary of knownGene cross-reference table by the table
    identifier, ie. refseq.
    """
    xref = {}
    for x in reader(kgxref):
        xref[x[table_id]] = x[table_symbol]
    return xref


def name_to_symbol(bed, kgxref, table_id, table_symbol):
    """Takes bed from a table, ie. refseq or knowngene, and converts the name
    field to more meaningful gene symbol.
    """
    xref = cross_ref(kgxref, table_id, table_symbol)
    bed = BedTool(bed)
    for b in bed:
        
        if "_" in b.chrom: continue
        
        # refseq beds have a bunch of descriptive garbage included on the name
        try:
            name = xref[b.name]
        except KeyError, e:
            sys.stderr.write(">> %s: Not found. Skipping...\n" % e)
        
        fields = (b.chrom, b.start, b.stop, name, b.score, b.strand)
        print "\t".join(map(str, fields))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('bed', help="input bed with gene id")
    p.add_argument('xref', help="cross-reference table from UCSC")
    p.add_argument('--id', dest="tableid", help="column in xref corresponding id")
    p.add_argument('--symbol', dest='tablesymbol', help="column in xref corresponding to symbol")
    args = p.parse_args()
    if not args.tableid or not args.tablesymbol:
        sys.exit(p.print_help())
    name_to_symbol(args.bed, args.xref, args.tableid, args.tablesymbol)


if __name__ == "__main__":
    sys.exit(main())