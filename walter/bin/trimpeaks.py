#!/usr/bin/env python
# encoding: utf-8
"""
Trim peaks by finding the point of maximum intensity then expanding to set
width.
"""
import os
import re
import sys
import subprocess as sp
import numpy as np
import tempfile
from genomedata import Genome
from pybedtools import BedTool
from toolshed import reader, header

def trim_peaks(results, genomedata, width, reference):
    # chr1:545523-545821
    ref = BedTool(reference)
    header = "chrom start stop baseMean baseMeanA baseMeanB foldChange log2FoldChange pval padj trimmedPeak sequence".split()
    print "\t".join(header)
    genome = Genome(genomedata)
    for toks in reader(results, sep=" "):
        if toks['padj'] == 'NA' or float(toks['padj']) > 0.01: continue
        chrom, coords = toks['id'].split(":")
        start, stop = coords.split("-")
        start = int(start)
        stop = int(stop)
        chromosome = genome[chrom]
        counts = chromosome[start:stop]
        summit = get_summit(start, counts)
        if not summit: continue
        sumstart = summit - width / 2
        sumstop = summit + width / 2
        sequence = chromosome.seq[sumstart:sumstop].tostring()
        fields = [chrom, start, stop, toks['baseMean'], toks['baseMeanA'],
                    toks['baseMeanB'], toks['foldChange'], toks['log2FoldChange'], 
                    toks['pval'], toks['padj'], "%s:%s-%s" % (chrom, sumstart, sumstop), 
                    sequence]
        print '\t'.join(map(str, fields))

def get_summit(start, counts):
    """Finds median position of the maximum intensity value."""
    maximum = 0
    maxpositions = []
    # get all of the positions with the max intensity
    for idx, count in enumerate(counts, start=start):
        # sum intensity across samples
        fullintensity = np.nansum(count)
        if fullintensity < maximum: continue
        if fullintensity == maximum:
            maxpositions.append(idx)
            continue
        if fullintensity > maximum:
            maximum = fullintensity
            del maxpositions[0:len(maxpositions)]
            maxpositions.append(idx)
    # get the middle position of the max counts
    if len(maxpositions) == 0:
        return None
    else:
        midposition = int(maxpositions[len(maxpositions) / 2])
        return midposition

def main(args):
    trim_peaks(args.results, args.genomedata, args.width, args.reference)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('results', help='DESeq results file')
    p.add_argument('genomedata', help="genome data archive")
    p.add_argument('reference', help="gtf")
    p.add_argument("-w", "--width", dest='width', default=60, 
                    type=int, help="full peak width [ 60 ]")
    args = p.parse_args()
    main(args)