#!/usr/bin/env python
# encoding: utf-8
"""
Add a piece or all of the index read back onto R1 (or R2).
"""
from itertools import izip, ifilter, islice, imap
from toolshed import nopen

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

def main(args):
    for idx, rec in izip(readfq(args.index), readfq(args.fastq)):
        assert idx.name.partition(" ")[0] == rec.name.partition(" ")[0]
        if args.end:
            idx.seq = idx.seq[args.begin:args.end]
            idx.qual = idx.qual[args.begin:args.end]
        else:
            idx.seq = idx.seq[args.begin:]
            idx.qual = idx.qual[args.begin:]
        rec.seq = idx.seq + rec.seq
        rec.qual = idx.qual + rec.qual
        print rec

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("index", metavar="INDEX", help="index read")
    p.add_argument("fastq", metavar="FASTQ", help="R1 or R2")
    p.add_argument("-b", dest="begin", default=0, type=int,
            help="0-based start of index read to save [%(default)s]")
    p.add_argument("-e", dest="end", type=int,
            help="0-based end of index read to save")
    args = p.parse_args()
    main(args)