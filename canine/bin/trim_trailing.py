#!/usr/bin/env python
# encoding: utf-8
"""
"""
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

def readfq(fq):
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield Fastq(rd)


def chunker(it, n):
    """
    >>> chunker("ABCDEF", 4)
    ['ABCD', 'BCDE', 'CDEF', 'DEF', 'EF', 'F']
    """
    return [it[i:i+n] for i in xrange(0, len(it), 1)]

def get_trim_loc(seq, base, fraction, wsize):
    for i, subseq in enumerate(chunker(seq, wsize)):
        if float(subseq.count(base))/len(subseq) >= fraction:
            return subseq.find(base) + i
    return len(seq)

def main(args):
    for rd in readfq(args.fastq):
        pos = get_trim_loc(rd.seq, args.b, args.f, args.w)
        if pos > args.l:
            rd.seq = seq[:pos]
            rd.qual = qual[:pos]
            print rd

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq")
    p.add_argument("-b", choices=['A','C','G','T','N'], required=True,
            help="base to trim")
    p.add_argument("-f", metavar="FRACTION", default=0.6, type=float,
            help="trimming threshold [%(default)s]")
    p.add_argument("-l", metavar="LENGTH", default=30, type=int,
            help="minimum sequence length [%(default)s]")
    p.add_argument("-w", metavar="WINDOW", default=15, type=int,
            help="window size [%(default)s]")
    args = p.parse_args()
    main(args)
