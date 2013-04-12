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
chr13:45,447,283-45,448,282 -- need 4 separate peaks here!
"""

__version__ = "0.1"

import os
import sys
import tempfile
import subprocess as sp
from toolshed import reader

def summit_start(starts, counts):
    """finds max position of starts given counts.
    
    >>> summit_start("20217423,20217424,20217425,20217426", "1,2,33,18")
    (20217425, 33)
    """
    starts = map(int, starts.split(","))
    counts = map(int, counts.split(","))
    summit_height = max(counts)
    return starts[counts.index(summit_height)], summit_height

def get_counts(bam):
    # tmp = tempfile.mkstemp(suffix=".bedgraph")[1]
    tmp = "tempcov.bedgraph"
    cmd = "bedtools genomecov -bg -5 -ibam %s > %s" % (bam, tmp)
    sp.call(cmd, shell=True)
    return tmp

def get_summits(bed, bedgraph):
    """
    peaks have no strand info, so it may be best to recount w/o respect to strand
    """
    # tmp = open(tempfile.mkstemp(suffix=".bed")[1], 'w')
    tmp = open("tmp_summits.bed", 'w')
    cmd = "|bedtools map -c 2 -o collapse -a %s -b %s |\
                 bedtools map -c 4 -o collapse -a - -b %s" % \
                 (bed, bedgraph, bedgraph)
    res = ["chr", "start", "stop", "name", "q", "starts", "counts"]
    for l in reader(cmd, header=res):
        # peak shifting isn't perfect so you'll get some peaks
        # where no reads overlap. those peaks typically have very few reads, so
        # just filter them out.
        if l['starts'] == ".":
            continue
        start, height = summit_start(l['starts'], l['counts'])
        fields = [l['chr'], start, start + 1, l['name'], l['q'], height]
        tmp.write("\t".join(map(str, fields)) + "\n")
    tmp.close()
    return tmp.name

def add_slop(bed, sizes, slop):
    """adds slop to bed entries. returns temp file name."""
    # tmp = tempfile.mkstemp(suffix=".bed")[1]
    tmp = "tempslop.bed"
    cmd = "bedtools slop -b %d -i %s -g %s > %s" % (slop, bed, sizes, tmp)
    sp.call(cmd, shell=True)
    return tmp

def classify_peaks(bed):
    """not sure yet"""
    # first need to add the slop to the summits bed
    # need the sequences for these regions
    pass

def main(args):
    print >>sys.stderr, ">> getting counts"
    counts_tmp = get_counts(args.bam)
    print >>sys.stderr, ">> locating summit per peak"
    summit_tmp = get_summits(args.bed, counts_tmp)
    print >>sys.stderr, ">> adding slop to summits"
    slop_tmp = add_slop(summit_tmp, args.sizes, args.wsize)
    print >>sys.stderr, ">> retrieving peak sequences"
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