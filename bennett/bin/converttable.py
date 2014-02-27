#!/usr/bin/env python
# coding=utf-8
"""
convert table to fasta and metadata table
"""

from argparse import ArgumentParser, RawDescriptionHelpFormatter


def main(table):
    fasta = open("seqs.fasta", 'w')
    metadata = open("meta.txt", 'w')
    print >>metadata, "#OTU ID\tsample"
    for line in open(table):
        toks = line.strip().split()
        # name
        print >>fasta, ">" + toks[0]
        # sequence
        print >>fasta, toks[1]
        # metadata output
        print >>metadata, "\t".join(toks)
    fasta.close()
    metadata.close()


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    p.add_argument("table")
    args = p.parse_args()
    main(args.table)
