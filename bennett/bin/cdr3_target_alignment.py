#!/usr/bin/env python
# encoding: utf-8
"""
Align observed CDR3 sequences to target sequences obtained outside of NGS.
"""

import sys
import argparse
from Bio import pairwise2
from toolshed import reader
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


class Alignment(object):
    __slots__ = ['query', 'target', 'score', 'start', 'stop']

    def __init__(self, args):
        for a, v in zip(self.__slots__, args):
            setattr(self, a, v)

    def __repr__(self):
        return "Alignment ({query} to {target})".format(query=self.query,
                target=self.target)

    def __str__(self):
        return "\n".join(str(getattr(self, s)) for s in self.__slots__)

    @property
    def alignment_length(self):
        return abs(self.stop - self.start)

def process_chunk(chunk):
    query, target = chunk
    print target
    return False
    for aln in pairwise2.align.localms(query, target, 1, -1, -1, -1):
        return Alignment(aln), target

def main(queries, targets, mismatches, shortest=5):
    queries = set([toks['CDR3-IMGT'] for toks in reader(queries) if toks['Functionality'] == "productive" and len(toks['CDR3-IMGT']) > shortest])
    targets = set([target.strip("\r\n") for target in open(targets, 'rb') if len(target) > shortest])

    print >>sys.stderr, "comparing", len(queries), "unique query sequences"
    print >>sys.stderr, "against", len(targets), "targets"

    for query in queries:
        for target in targets:
            for a_query, a_target, score, start, stop in pairwise2.align.localms(query, target, 1, -1, -5, -1):
                if score > len(target) - mismatches:
                    print a_query
                    print a_target

if __name__ == '__main__':
    p = ArgumentParser(description=__doc__,
            formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('queries', help="IMGT AA sequence file (5)")
    p.add_argument('targets', help="target sequences")
    p.add_argument("-m", "--mismatches", type=int, default=3)
    args = vars(p.parse_args())

    main(**args)
