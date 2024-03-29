#!/usr/bin/env python
# encoding: utf-8
"""
Add a piece or all of the index read back onto R1.
"""

from itertools import izip, islice
from toolshed import nopen


def readfq(fq):
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


    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield Fastq(rd)


def main(index, fastq, begin, length):
    end = begin + length
    for idx, rec in izip(readfq(index), readfq(fastq)):
        assert idx.name.partition(" ")[0] == rec.name.partition(" ")[0]
        idx.seq = idx.seq[begin:end]
        idx.qual = idx.qual[begin:end]
        assert len(idx.seq) == length
        rec.seq = idx.seq + rec.seq
        rec.qual = idx.qual + rec.qual
        print rec


if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("index", metavar="INDEX", help="index read")
    p.add_argument("fastq", metavar="FASTQ", help="R1 fastq")
    p.add_argument("-b", dest="begin", default=0, type=int,
            help="inclusive 0-based start of UMI")
    p.add_argument("-l", dest="length", type=int, default=None,
            help="length of UMI")
    args = vars(p.parse_args())
    main(**args)
