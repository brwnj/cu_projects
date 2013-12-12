#!/usr/bin/env python
# encoding: utf-8
"""
From the positive and negative peak files, calculate the average number of 
feature peaks per gene.
"""
import argparse
import sys
from pybedtools import BedTool
from toolshed import reader
from collections import defaultdict, Counter

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def get_peak_counter(bed, ref_bed):
    """
    bed is either the pos or neg strand file of peak coordinates
    ref_bed is the reference annotation for exon, intron, etc.
    
    returns defaultdict(Counter)
    """
    bed = BedTool(bed)
    ref_bed = BedTool(ref_bed)
    peaks = defaultdict(Counter)
    for peak in bed.intersect(ref_bed,f=0.5,wo=True):
        gene_name = peak[6]
        peaks[gene_name].update([gene_name])
    return peaks


def get_total(defdict, filename):
    """from defaultdict(Counter)s, return the sum of all counters. if a
    filename is provided, write the values to a file.
    
    returns float
    """
    if filename:
        fout = open('%s.txt' % filename, 'w')
    
    total = 0.0
    for gene in defdict.keys():
        assert len(defdict[gene].values()) == 1
        if filename:
            fout.write("%s\n" % defdict[gene].values()[0])
        total += defdict[gene].values()[0]
    return total


def main():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('pos', help='peaks on the positive strand')
    p.add_argument('neg', help='peaks on the negative strand')
    p.add_argument('ref', help='annotation file')
    p.add_argument('name', help='name of annotation feature')
    p.add_argument('-o', '--output', help='file name prefix for counts')
    args = p.parse_args()
    
    positive = get_peak_counter(args.pos, args.ref)
    negative = get_peak_counter(args.neg, args.ref)

    total = get_total(positive, args.output) + \
                    get_total(negative, args.output)

    if args.output:
        total = get_total(positive, args.output) + \
                    get_total(negative, args.output)
    else:
        total = get_total(positive, False) + get_total(negative, False)
    
    print "Peaks per mRNA (feature = %s): %.0f" % \
                (args.name, (total / (len(positive) + len(negative))))


if __name__ == "__main__":
    main()