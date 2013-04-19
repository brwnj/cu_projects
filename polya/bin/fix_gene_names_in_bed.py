#!/usr/bin/env python
# encoding: utf-8
"""
Give the same gene name across the same polya "family".
"""
from itertools import groupby
from toolshed import nopen

def read_polya_bed(fh):
    for key, grp in groupby(fh, key=\
                lambda item: item.split("\t")[3].split(":")[1].split(".")[1]):
        yield grp

def main(args):
    head = ["chrom", "start", "stop", "name", "score", "strand"]
    with nopen(args.bed) as bed:
        for group in read_polya_bed(bed):
            gene = None
            for line in group:
                d = dict(zip(head, line.strip().split()))
                if gene is None:
                    gene = d['name'].split(":")[0]
                d['name'] = "%s:%s" % (gene, d['name'].split(":")[1])
                print "\t".join([d[h] for h in head])

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bed",
            help="bed file with different gene names across polya sites.")
    main(p.parse_args())
