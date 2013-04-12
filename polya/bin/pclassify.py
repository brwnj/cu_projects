#!/usr/bin/env python
# encoding: utf-8
"""
have peaks bed
data in the genomedata archive is coming from bedgraphs where counts were

for each summit; grab 100 nt upstream

category
1   has A[A,T]TAAA; NOT A-rich downstream from this cleavage site
2   has A[A,T]TAAA; with A-rich sequence downstream
3   lacks A[A,T]TAAA; no A-rich region downstream
4   only downstream A-rich sequence

100 bp from known polya sites

filter out peaks with less than 10 reads of support
canonical PAS should be found in -10 to -30 range of the cleavage site
no more than three canonical PAS should be located in the upstream window
chr6:90,095,714-90,096,087
"""

__version__ = "0.1"

import os
import sys
import tempfile
import subprocess as sp
from toolshed import nopen

def get_summit(bed, bam, sizes):
    """
    peaks have no strand info, so it may be best to recount w/o respect to strand
    """
    print >>sys.stderr, "getting counts"
    # tmp = tempfile.mkstemp(suffix=".count")[1]
    tmp = "tempsummit.bed"
    cmd = "bedtools genomecov -bg -5 -ibam %s -g %s | \
                bedtools map -a %s -b - -c 4 -o max > %s" % \
                (bam, sizes, bed, tmp)
    sp.call(cmd, shell=True)
    return tmp

def add_slop(bed, sizes, slop):
    """adds slop to bed entries. returns temp file name."""
    # tmp = open(tempfile.mkstemp(suffix=".bed")[1], 'w')
    print >>sys.stderr, "adding slop to summit"
    tmp = open("tempslop.bed", 'w')
    cmd = "|bedtools slop -b %d -i %s -g %s" % (slop, bed, sizes)
    tmp.write("".join([line for line in nopen(cmd)]))
    tmp.close()
    return tmp.name

def classify_peaks(bed):
    """bed file represents summit locations."""
    # first need to add the slop to the summits bed
    # need the sequences for these regions
    pass

def main(args):
    # summits is bullshit
    # find the highest point of the peak; aka THE SUMMIT!
    # summit_tmp = get_summit(args.bed, args.bam, args.sizes)
    summit_tmp = "tempsummit.bed"
    # bed to larger regions
    slop_tmp = add_slop(args.bed, args.sizes, args.wsize)
    # os.remove(slop_tmp)
    # larger regions to fasta
    # classify the lines of the fasta
    # classify_peaks()

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__, version=__version__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('bed', metavar="BED", help="bed file of peak regions")
    p.add_argument('bam', metavar="BAM", help="aligned reads")
    p.add_argument('fasta', metavar="FASTA", help="genomic reference sequence")
    p.add_argument('sizes', metavar="SIZES", help="genomic chromosome sizes")
    p.add_argument('-w', '--window-size', dest="wsize", type=int, default=50,
            help="number of nucleotides to add up- and downstream of PAS [%(default)s]")
    main(p.parse_args())