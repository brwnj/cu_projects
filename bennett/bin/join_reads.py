#!/usr/bin/env python
# encoding: utf-8
"""
Join reads based on local alignment, taking higher quality base where mismatches
are present.
"""
import sys, string
from Bio import pairwise2
from toolshed import nopen
from itertools import islice, izip

COMPLEMENT = string.maketrans('ACGTNSRYMKWHBVD','TGCANSRYMKWHBVD')

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

def rev_comp(seq):
    """return reverse complement of seq."""
    return seq.translate(COMPLEMENT)[::-1]

def decode(x):
    return ord(x) - 33

def main(args):
    """For each pair, perform local alignment, fix the mismatches, and output
    joined read with qual.
    """
    pw = pairwise2.align.localms
    for (r1, r2) in izip(readfq(args.R1), readfq(args.R2)):
        r1.name = r1.name.split()[0]
        r2.name = r2.name.split()[0]
        assert r1.name == r2.name
        # reverse r2
        r2.seq = rev_comp(r2.seq)
        r2.qual = r2.qual[::-1]
        try:
            # get best match
            r1aln, r2aln, score, start, stop = pw(r1.seq, r2.seq, 1, -1, -5, -3)[0]
            # assume overlapping is NOT occurring on 5' end
            if r1aln.startswith("-----"): continue
            seq = ""
            qual = ""
            idx1 = 0
            idx2 = 0
            for (a, b) in izip(r1aln, r2aln):
                aq = r1.qual[idx1]
                bq = r2.qual[idx2]
                if a == "-":
                    seq += b
                    qual += bq
                    idx2 += 1
                elif b == "-":
                    seq += a
                    qual += aq
                    idx1 += 1
                elif a == b:
                    seq += a
                    if decode(aq) > decode(bq):
                        qual += aq
                    else:
                        qual += bq
                else:   # a != b and they're both bases
                    if decode(aq) > decode(bq):
                        seq += a
                        qual += aq
                    else:
                        seq += b
                        qual += bq
            assert len(seq) == len(qual)
            r1.seq = seq
            r1.qual = qual            
            print r1
        except IndexError:  # no overlap exists
            continue

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("R1", help="fastq")
    p.add_argument("R2", help="fastq")
    args = p.parse_args()
    main(args)
