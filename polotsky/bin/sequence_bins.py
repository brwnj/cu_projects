#!/usr/bin/env python
# encoding: utf-8
"""
Finds and quantifies unique and similar sequences among a FASTQ. Entire FASTQ
will be read into memory.
"""
import sys
import editdist as ed
from collections import Counter
from toolshed import nopen
from itertools import islice

def read_fastq(fh):
    """fastq parser that returns name, seq, and qual."""
    while True:
        values = list(islice(fh, 4))
        if len(values) == 4:
            id1, seq, id2, qual = values
        elif len(values) == 0:
            raise StopIteration
        else:
            raise EOFError("unexpected end of file")
        assert id1.startswith('@')
        assert id2.startswith('+')
        assert len(seq) == len(qual)
        yield id1[1:-1], seq[:-1], qual[:-1]

def match(target, query, mismatches):
    """Run local alignment. Returns bool."""
    try:
        best_score = pairwise2.align.localms(target, query, 1, -1, -2, -1)[0][2]
        if best_score >= (len(query) - mismatches):
            return True
    except IndexError:
        return False

def group_unique(fastq):
    """Group identical reads using a Counter. Returns Counter."""
    c = Counter()
    with nopen(fastq) as fh:
        for name, seq, qual in read_fastq(fh):
            c.update([seq])
    return c

def group_matches(counter, mismatches):
    """not sure yet..."""
    d = {}
    seqs = list(counter).sort(key = len)
    print seqs
    # ignore = []
    #     
    # order the sequences by length
    # iterate over the ordered sequences
    # find anything that is a match and add to ignore

    # ignore = set()
    # seqs = set(list(counter))
    # for query in seqs:
    #     if query in ignore: continue
    #     for target in seqs:
    #         if target == query: continue
    #         if ed.distance(query, target) < mismatches:
    #             ignore.add(target)
    # return seqs - ignore

def main(args):
    print >>sys.stderr, ">> grouping unique reads"
    reads = group_unique(args.fastq)
    if args.mismatches == 0:
        for seq, count in reads.iteritems():
            print "%s\t%d" % (seq, count)
    print >>sys.stderr, ">> further grouping allowing mismatches"
    reads = group_matches(reads, args.mismatches)
    

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq", metavar="FASTQ", help="reads to process")
    p.add_argument("-m", "--mismatches", type=int, default=3,
            help="number of mismatches allowed when combining bins [%(default)s]")
    main(p.parse_args())