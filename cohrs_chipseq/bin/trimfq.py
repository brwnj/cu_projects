#!/usr/bin/env python
# encoding: utf-8
"""
trim a fastq in various ways. execution of options occurs in the order they are
listed.
"""

import argparse
from toolshed import nopen
from itertools import islice

class Fastq(object):
    def __init__(self, args):
        self.name = args[0][1:]
        self.seq = args[1]
        self.qual = args[3]
        assert len(self.seq) == len(self.qual)

    def __repr__(self):
        return "Fastq({name})".format(name=self.name)

    def __str__(self):
        return "@{name}\n{seq}\n+\n{qual}".format(name=self.name,
                seq=self.seq, qual=self.qual)

    def trim_max(self, n):
        self.seq = self.seq[:n]
        self.qual = self.qual[:n]
    
    def trim_left(self, n):
        self.seq = self.seq[n:]
        self.qual = self.qual[n:]
    
    def trim_right(self, n):
        self.seq = self.seq[:-n]
        self.qual = self.qual[:-n]

def readfq(fq):
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield Fastq(rd)

def main(fastq, m, b, e, length_threshold):
    fq = readfq(fastq)
    for f in readfq(fastq):
        if m != 0:
            f.trim_max(m)
        if b != 0:
            f.trim_left(b)
        if e != 0:
            f.trim_right(e)
        if len(f.seq) >= length_threshold:
            print f

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('fastq')
    p.add_argument('-m', metavar="INT", type=int, default=0, help="trim from right to maximum of INT bp")
    p.add_argument('-b', metavar="INT", type=int, default=0, help="trim INT bp from left")
    p.add_argument('-e', metavar="INT", type=int, default=0, help="trim INT bp from right")
    p.add_argument('-l', metavar="INT", type=int, default=18, help="shortest allowable read length")
    args = p.parse_args()
    main(args.fastq, args.m, args.b, args.e, args.l)
