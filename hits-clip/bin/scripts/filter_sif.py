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
    for a in reader(attrs_file, header="mir op value".split(), sep=' '):
        try:
            mir_filter[a['mir']] = a['value']
        except KeyError:
            pass
    return mir_filter


def main():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('sif', help='unfiltered network (sif)')
    p.add_argument('attrs', help='Cytoscape attributes file for abundance')
    args = p.parse_args()
    
    mirna_filter = get_dict(args.attrs)
    
    for a in reader(args.sif, header="mir length gene".split(), sep='\t'):
        try:
            m_filter = mirna_filter[a['mir']] is "0.0"
        except KeyError:
            # header or miRNA present in the mature.fa but not the gff3, ignore for now
            continue
        if m_filter: continue
        fields = (a['mir'], a['length'], a['gene'])
        print "\t".join(map(str, fields))


if __name__ == "__main__":
    main()