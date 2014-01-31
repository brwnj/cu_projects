#!/usr/bin/env python
# encoding: utf-8
"""
CCGCTGGAAGTGACTGACAC

rev comp
GTGTCAGTCACTTCCAGCGG
"""

import sys
from numpy import average
from toolshed import nopen
from itertools import islice
from string import maketrans
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


COMPLEMENT = maketrans('ATCG', 'TAGC')


def revcomp(seq):
    return seq.translate(COMPLEMENT)[::-1]


def readfq(fq):
    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            r = [x for x in islice(fqclean, 4)]
            if not r: raise StopIteration
            assert all(r) and len(r) == 4
            # seq and qual length
            assert len(r[1]) == len(r[3])
            yield r[0], r[1], r[3]


def print_read(name, seq, qual):
    print "{name}\n{seq}\n+\n{qual}".format(name=name, seq=seq, qual=qual)


def main(fastq, adapter_sequence, seed_length, min_read_length, reverse_complement):
    if reverse_complement:
        adapter_sequence = revcomp(adapter_sequence)[:seed_length]
    else:
        adapter_sequence = adapter_sequence[:seed_length]
    total_reads = 0
    index_start_locations_passed = []
    index_start_locations_failed = []
    seed_not_found = 0
    for name, seq, qual in readfq(fastq):
        total_reads += 1

        try:
            adapter_start = seq.index(adapter_sequence)
            if adapter_start > min_read_length:
                index_start_locations_passed.append(adapter_start)
                print_read(name, seq[:adapter_start], qual[:adapter_start])
            else:
                index_start_locations_failed.append(adapter_start)
        except ValueError:
            # full seed not found
            # see if the read ends with at least 3 bases of adapter sequence
            seed_found = False

            for i in xrange(seed_length-1, 2, -1):
                if seq.endswith(adapter_sequence[:i]):
                    seed_found = True
                    adapter_start = len(seq) - i
                    index_start_locations_passed.append(adapter_start)
                    print_read(name, seq[:adapter_start], qual[:adapter_start])

            if not seed_found:
                seed_not_found += 1

    print >>sys.stderr, fastq
    print >>sys.stderr, "Total reads:", total_reads
    print >>sys.stderr, "Passing reads:", len(index_start_locations_passed)
    print >>sys.stderr, "Average passing read length:", average(index_start_locations_passed)
    print >>sys.stderr, "Failed reads due to length:", len(index_start_locations_failed)
    print >>sys.stderr, "Average failed read length:", average(index_start_locations_failed)
    print >>sys.stderr, "Failed reads due to missing seed sequence:", seed_not_found


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('fastq')
    p.add_argument('-a', '--adapter-sequence')
    p.add_argument('-l', '--seed-length', type=int, default=6, help="bases of adapter to match")
    p.add_argument('-m', '--min-read-length', type=int, default=18, help="minimum allowable read length after trimming")
    p.add_argument('--reverse-complement', action="store_true",
            help="acts upon the adapter sequence")
    args = vars(p.parse_args())
    main(**args)
