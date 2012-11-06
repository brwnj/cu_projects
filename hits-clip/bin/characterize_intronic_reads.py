#!/usr/bin/env python
# encoding: utf-8
"""
Find the distance each intronic read falls from either the 3' or 5' UTR.
"""
import argparse
import sys
from pybedtools import BedTool

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def main():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('reads', help='bed format file of reads')
    p.add_argument('intron', help='intron coordinates bed')
    p.add_argument('utr', help="UTR coordinates bed")
    p.add_argument('-v', '--verbose', action='store_true',
            help="maximum verbosity")

    args = p.parse_args()
    
    bed = BedTool(args.reads)
    intron = BedTool(args.intron)
    utr = BedTool(args.utr)
    
    annotated_bed = {}
    
    if args.verbose: print >>sys.stderr, ">> finding intersection"
    intersection = bed.intersect(intron, s=True)
    
    for row in intersection:
        annotated_bed[row.name] = {'chrom':row.chrom, 'start':row.start, \
                                        'stop':row.stop}
    
    if args.verbose: print >>sys.stderr, ">> annotating"
    
    for row in intersection.closest(utr, s=True, d=True):
        annotated_bed[row.name]['utr'] = int(row.fields[-1:][0])

    step = 100
    rangemax = 5000
    
    bin_counts = {}
    max_utr = 0
    
    for ubound in range(step, rangemax + 1, step):
        bincount = 0
        for name, val in annotated_bed.iteritems():
            if val['utr'] <= ubound and val['utr'] > (ubound - step):
                bincount += 1
            if val['utr'] > max_utr:
                max_utr = val['utr']
        fields = (ubound, bincount)
        print "\t".join(map(str, fields))
    
    for ubound in range(50000, max_utr, 50000):
        bincount = 0
        for name, val in annotated_bed.iteritems():
            if val['utr'] <= ubound and val['utr'] > (ubound - 100000):
                bincount += 1
        fields = (ubound, bincount)
        print "\t".join(map(str, fields))


if __name__ == "__main__":
    main()