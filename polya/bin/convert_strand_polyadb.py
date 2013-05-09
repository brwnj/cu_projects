#!/usr/bin/env python
# encoding: utf-8
"""
convert strand of polya_db to their actual strand.
"""

from toolshed import nopen
from itertools import groupby

def polya_site(item):
    return item.split(".")[2]

def polya_name(item):
    return item.split("\t")[3].split(".")[1]

def read_polya_bed(fh):
    for key, grp in groupby(fh, key=lambda item: polya_name(item)):
        yield grp

def get_strand(sites):
    i = None
    if len(sites) == 1:
        return "+"
    for site in sites:
        if i is None:
            i = site
            continue
        return "+" if site > i else "-"

def main(args):
    head = ["chrom", "start", "stop", "name", "score", "strand"]
    with nopen(args.bed) as bed:
        for group in read_polya_bed(bed):
            sites = []
            grp = []
            for l in group:
                l = dict(zip(head, l.strip().split()))
                grp.append(l)
                sites.append(polya_site(l['name']))
            strand = get_strand(sites)
            for l in grp:
                l['strand'] = strand
                print "\t".join([l[h] for h in head])

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bed", help="bed of polyadb from ucsc")
    args = p.parse_args()
    main(args)
