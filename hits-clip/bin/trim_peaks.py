#!/usr/bin/env python
# encoding: utf-8
"""
Trim peaks by finding the point of maximum intensity then expanding to the
size of width.
"""
import argparse
import os
import sys
from genomedata import Genome
from pybedtools import BedTool
from numpy import nansum


def trim_peaks(bed, genomedata, track, chrom, width):
    """For each peak, find mid then write new coordinates from the mid.
    """
    bed = BedTool(bed)
    
    with Genome(genomedata) as genome:
        for b in bed:
            # parallelize by chrom if supplied
            if chrom and b.chrom != chrom: continue
            
            chromosome = genome[b.chrom]
            
            # intensities across specified tracks for peak region
            intensities = chromosome[b.start:b.stop, track]
            mid = get_mid(b.chrom, b.start, intensities)
            
            # calling intensities for some regions were returning all 'nan'
            if not mid: continue
            
            start = mid - width / 2
            stop = mid + width / 2
            
            # output new peak coordinates
            fields = (b.chrom, start, stop)
            print '\t'.join(map(str, fields))


def get_mid(chrom, start, intensities):
    """Finds median position of the maximum intensity value.
    """
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
    
    # get the middle position of the max intensities
    if len(maxpositions) == 0:
        return None
    else:
        midposition = int(maxpositions[len(maxpositions) / 2])
        return midposition


def main():
    p = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('peaks', help='peaks file (bed)')
    p.add_argument('genomedata', help="genome data archive")
    p.add_argument('trackname', nargs='+', 
                    help='space delimited track names')
    p.add_argument("-w", "--width", dest='width', default=60, 
                    help="full peak width", type=int)
    p.add_argument("-c", "--chrom", dest='chromosome', default=None, 
                    help="parallelize by specifying chromosome name")
    
    args = p.parse_args()
    trim_peaks(args.peaks, args.genomedata, args.trackname, \
                    args.chromosome, args.width)


if __name__ == "__main__":
    sys.exit(main())