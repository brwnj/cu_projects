#!/usr/bin/env python
# encoding: utf-8
"""
Finds and quantifies unique and similar sequences among a FASTQ. Entire FASTQ
will be read into memory.
"""
import sys
import editdist as ed
from toolshed import nopen
from itertools import islice
from collections import Counter

def read_fastq(fh):
    """FASTQ parser that yields name, seq, and qual."""
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

def group_unique(fastq):
    """Group identical reads using a Counter. Returns Counter."""
    c = Counter()
    with nopen(fastq) as fh:
        for name, seq, qual in read_fastq(fh):
            c.update([seq])
    return c

def distance(a, b):
    """Calculate Levenshtein distance accounting for length differences between
    the two strings. Returns int.
    >>> distance('CATGGGTGGTTCAGTGGTAGAATTCTCGCCTGCC', 'GTGCTGTAGGCATT')
    2
    """
    return ed.distance(a, b) - abs(len(a) - len(b))

def group_matches(counter, mismatches):
    """not sure yet..."""
    seen = set()
    # ordered by length to return longest sequence
    seqs = list(counter)
    seqs.sort(key = len, reverse = True)
    for target in seqs:
        if target in seen: continue
        seen.update([target])
        for query in seqs:
            if query in seen: continue
            # check mismatch tolerance
            if distance(target, query) < mismatches:
                counter[target] += counter[query]
                # set added items to zero to mark for removal
                counter[query] = 0
                # and don't compare against this sequence any more
                seen.update([query])
    # remove 0s from the counter
    counter += Counter()
    return counter

def main(args):
    print >>sys.stderr, ">> grouping unique reads"
    reads = group_unique(args.fastq)
    if args.mismatches == 0:
        for seq, count in reads.iteritems():
            print "%s\t%d" % (seq, count)
    else:
        print >>sys.stderr, \
                ">> accounting for %d mismatch tolerance" % args.mismatches
        reads = group_matches(reads, args.mismatches)
        for seq, count in reads.iteritems():
            print "%s\t%d" % (seq, count)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq", metavar="FASTQ", help="reads to process")
    p.add_argument("-m", "--mismatches", type=int, default=3,
            help="number of mismatches allowed when combining bins [%(default)s]")
    main(p.parse_args())