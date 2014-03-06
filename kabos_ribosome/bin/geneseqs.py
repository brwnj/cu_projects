#!/usr/bin/env python
# coding=utf-8
"""
"""

import sys
from itertools import groupby
from toolshed import nopen, reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def readfa(fa):
    class Fasta(object):
        def __init__(self, name, seq):
            self.name = name
            self.seq = seq

        def __repr__(self):
            return "Fasta({name})".format(name=self.name)

        def __str__(self, wrap=100):
            return ">{name}\n{seq}".format(name=self.name,
                                           seq="\n".join([self.seq[i:i + wrap] for i in range(0, len(self.seq), wrap)]))

    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip().split("|")[0]
            else:
                seq = ''.join(line.strip() for line in group)
                yield Fasta(name, seq)


def main(fasta, genelist):
    genes = set()
    for t in reader(genelist):
        genes.add(t['Ensembl Gene ID'])
    seen = set()
    for f in readfa(fasta):
        if f.name in genes:
            seen.add(f.name)
            print f
    print >>sys.stderr, "No sequence found for:"
    for g in genes - seen:
        print >>sys.stderr, g


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('fasta')
    p.add_argument('genelist')
    args = vars(p.parse_args())
    main(**args)
