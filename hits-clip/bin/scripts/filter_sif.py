#!/usr/bin/env python
# encoding: utf-8
"""
Filters sif to only contain miRNA interactions where the miRNA was shown to be
present.
"""
import argparse
import sys
from toolshed import reader

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def get_dict(attrs_file):
    """Return dictionary of miRNA abundance."""
    mir_filter = {}
    for toks in reader(attrs_file, header=["mir", "op", "value"], sep=' '):
        mir_filter[toks['mir']] = float(toks['value'])
    return mir_filter


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sif', help='unfiltered network (.sif)')
    p.add_argument('attrs', help='Cytoscape attributes file for abundance')
    args = p.parse_args()

    mirna_filter = get_dict(args.attrs)

    for a in reader(args.sif, header="mir length gene".split(), sep='\t'):
        try:
            m_filter = mirna_filter[a['mir']] == 0
        except KeyError:
            # header or miRNA present in the mature.fa but not the gff3, ignore for now
            continue
        try:
            if mirna_filter[a['mir']] == 0: continue
        except KeyError:
            continue
        else:
            print "\t".join([a['mir'], a['length'], a['gene']])


if __name__ == "__main__":
    main()
