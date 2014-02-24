#!/usr/bin/env python
# encoding: utf-8
"""
First: attempt to trim adapter sequence from read; if long enough after
identifying trim position, print read to stdout. Second: where no adapter was
found, search for discard search. When found, discard read. Third: if neither
were found, print read with altered read name (name:inspect).
"""

import sys
from numpy import average
from toolshed import nopen
from itertools import islice
from string import maketrans
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


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


def main(fastq, adapter_sequence, discard_sequence, min_read_length, mark):
    total_reads = 0
    index_start_locations_passed = []
    index_start_locations_failed = []
    seq_not_found = 0
    discard_found = 0

    for name, seq, qual in readfq(fastq):
        total_reads += 1

        # find the adapter sequence
        adapter_start = seq.find(adapter_sequence)
        if adapter_start > -1:
            if adapter_start > min_read_length:
                index_start_locations_passed.append(adapter_start)
                print_read(name, seq[:adapter_start], qual[:adapter_start])
            else:
                index_start_locations_failed.append(adapter_start)

        # find the discard sequence
        elif seq.find(discard_sequence) > -1:
            discard_found += 1

        # neither sequence found
        else:
            # if the read ends with either primer they will either be soft-
            # clipped or discarded in subsequent steps.

            # alter read name so we can later find this in the bam
            name += ":" + mark
            seq_not_found += 1
            print_read(name, seq, qual)

    print >>sys.stderr, fastq
    print >>sys.stderr, "Total reads:", total_reads
    passed_length = len(index_start_locations_passed)
    print >>sys.stderr, "Passing reads:", passed_length
    if passed_length > 0:
        print >>sys.stderr, "Average passing read length:", average(index_start_locations_passed)
    failed_length = len(index_start_locations_failed)
    print >>sys.stderr, "Failed reads due to length:", failed_length
    if failed_length > 0:
        print >>sys.stderr, "Average failed read length:", average(index_start_locations_failed)
    print >>sys.stderr, "Failed reads due to discard seq:", discard_found
    print >>sys.stderr, "Reads missing adapter sequence:", seq_not_found


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('fastq')
    p.add_argument('-a', '--adapter-sequence', required=True, help="adapter sequence to search for and trim read")
    p.add_argument('-d', '--discard-sequence', required=True, help="sequence to search for and discard entire read")
    p.add_argument('-m', '--min-read-length', type=int, default=16, help="minimum allowable read length after trimming")
    p.add_argument('-M', '--mark', default="inspect", help="string to add onto name; will later search for this among aligned reads")
    args = vars(p.parse_args())
    main(**args)
