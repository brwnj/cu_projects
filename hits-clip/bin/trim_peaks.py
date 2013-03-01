#!/usr/bin/env python
# encoding: utf-8
"""
Trim peaks by finding the point of maximum intensity then expanding to the
size of width.
"""

import os
import sys
from genomedata import Genome
from pybedtools import BedTool
from numpy import nansum

def trim_peaks(bed, genomedata, track, chrom, width, verbose):
    """For each peak, find mid then write new coordinates from the mid."""
    bed = BedTool(bed)

    # ideally this would look more like
    # for chromosome in genome:
    #     for supercontig, continuous in chromosome.itercontinuous():
    # but i've noticed data not being pulled out accurately
    
    with Genome(genomedata) as genome:
        for b in bed:
            # parallelize by chrom if supplied
            if chrom and b.chrom != chrom: continue
            # data accessor from genomedata
            chromosome = genome[b.chrom]
            # intensities across specified tracks for peak region
            intensities = chromosome[b.start:b.stop, track]
            # print intensities
            mid = get_mid(b.chrom, b.start, intensities)

            # if not mid: continue
            start = mid - width / 2
            stop = mid + width / 2
            # output new peak coordinates
            fields = (b.chrom, start, stop, b.name, b.score, b.strand)
            print '\t'.join(map(str, fields))

def get_mid(chrom, start, intensities):
    """Finds median position of the maximum intensity value."""
    maximum = 0
    maxpositions = []
    
    # get all of the positions with the max intensity
    for idx, intensity in enumerate(intensities, start=start):
        
        # sum intensity across samples
        fullintensity = nansum(intensity)
        if fullintensity < maximum: continue
        if fullintensity == maximum:
            maxpositions.append(idx)
            continue
        if fullintensity > maximum:
            maximum = fullintensity
            del maxpositions[0:len(maxpositions)]
            maxpositions.append(idx)
    
    # return midpoint among the max intensities
    if len(maxpositions) == 0:
        return None
    else:
        midposition = int(maxpositions[len(maxpositions) / 2])
        return midposition

def main(args):
    kwargs = {'bed':args.peaks,
              'genomedata':args.genomedata,
              'track':args.trackname,
              'chrom':args.chrom,
              'width':args.width,
              'verbose':args.verbose}
    trim_peaks(**kwargs)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('peaks', help='peaks file (bed)')
    p.add_argument('genomedata', help="genome data archive")
    p.add_argument('trackname', nargs='+', 
                    help='space delimited track names')
    p.add_argument("-w", "--width", default=60, type=int,
                    help="full peak width")
    p.add_argument("-c", "--chrom", default=None, 
                    help="parallelize by specifying chromosome name")
    p.add_argument("-v", "--verbose", action='store_true',
                    help="maximum verbosity")
    main(p.parse_args())