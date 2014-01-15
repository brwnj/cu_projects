#!/usr/bin/env python
# encoding: utf-8
"""
Align observed CDR3 sequences to target sequences obtained outside of NGS.
"""
import sys
import argparse
import multiprocessing
from itertools import product
from toolshed import reader
from Bio import pairwise2

# https://gist.github.com/aljungberg/626518
from multiprocessing.pool import IMapIterator

def wrapper(func):
  def wrap(self, timeout=None):
    # Note: the timeout of 1 googol seconds introduces a rather subtle
    # bug for Python scripts intended to run many times the age of the universe.
    return func(self, timeout=timeout if timeout is not None else 1e100)
  return wrap
IMapIterator.next = wrapper(IMapIterator.next)

class Alignment(object):
    __slots__ = ['query', 'target', 'score', 'start', 'stop']

    def __init__(self, args):
        for a, v in zip(self.__slots__, args):
            setattr(self, a, v)
        self.score = float(self.score)
        self.start = int(self.start)
        self.stop = int(self.stop)

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
    for aln in pairwise2.align.localms(query, target, 1, -1, -1, -1):
        return Alignment(aln), target

def main(cdr3s, targets, mismatches, threads):
    cdr3s = [cdr3.strip("\r\n") for cdr3 in open(cdr3s, 'rb')]
    targets = [target.strip("\r\n") for target in open(targets, 'rb')]

    p = multiprocessing.Pool(threads)
    for result in p.imap(process_chunk, product(cdr3s, targets), 100):
        try:
            alignment, target_seq = result
        except TypeError:
            # no alignments
            continue
        try:
            if alignment.score > len(target_seq) - mismatches:
                print alignment
        except AttributeError:
            # no alignments
            pass

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('observed_cdr3')
    p.add_argument('targets')
    p.add_argument("-m", "--mismatches", type=int, default=3)
    p.add_argument("-t", "--threads", type=int, default=1)
    args = vars(p.parse_args())

    cpus = multiprocessing.cpu_count()
    args['threads'] = cpus if args['threads'] > cpus

    main(**args)
