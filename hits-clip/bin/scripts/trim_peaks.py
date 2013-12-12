#!/usr/bin/env python
# encoding: utf-8
"""
Trim peaks by finding the point of maximum intensity then expanding to the
size of width.

Determines which tracks will be used for trimming based on file name and has 
assumptions about your genomedata archive:

1. your genomedata archive only has this sample once
2. sample names terminate with a period character, then list strand, e.g. 
MP1.pos.rmd
<sample>.<strand>.<...>

The assumptions concerning your peaks file is that it starts with the sample name
and is delimited by periods:

MP1.peaks.rmd.bed
<sample>.peaks.<...>
"""

import os
import sys
import warnings
from genomedata import Genome
from pybedtools import BedTool
from numpy import nansum

def trim_peaks(bed, genomedata, track_idx, peak_ext, verbose):
    """For each peak, find mid then write new coordinates from the mid."""
    bed = BedTool(bed)
    with Genome(genomedata) as gd:
        for chrom in gd:
            if verbose: sys.stderr.write(">> processing %s...\n" % chrom)
            for supercontig, continuous in chrom.itercontinuous():
                whole_tracks = {}
                # track_idx = {"+":1,"-":0}
                for strand, idx in track_idx.iteritems():
                    whole_tracks[strand] = continuous[:, idx]
                # will be iterating over this for each supercontig
                for b in bed:
                    # speed up processing
                    if b.chrom != chrom.name: continue
                    # only test the correct peaks
                    # should work even if not sorted
                    if not b.start >= supercontig.start: continue
                    if not b.stop <= supercontig.end: continue
                    # get the relative position since continous will start at 0
                    relative_start = supercontig.project(b.start)
                    relative_stop = supercontig.project(b.stop)
                    # counts are only taken from appropriate strand
                    counts = whole_tracks[b.strand][relative_start:relative_stop]
                    # pass the actual start since we're trimming the peak here
                    mid = get_mid(b.start, counts)
                    start = mid - peak_ext
                    stop = mid + peak_ext
                    fields = [b.chrom, start, stop, b.name, b.score, b.strand]
                    print "\t".join(map(str, fields))

def get_mid(start, counts):
    """Finds median position of the maximum intensity value."""
    maximum = 0
    maxpositions = []
    # get all of the positions with the max intensity
    for idx, count in enumerate(counts, start=start):
        if count < maximum: continue
        if count == maximum:
            maxpositions.append(idx)
            continue
        if count > maximum:
            maximum = count
            del maxpositions[0:len(maxpositions)]
            maxpositions.append(idx)
    # return midpoint among the max intensities
    return int(maxpositions[len(maxpositions) / 2])

def get_track_idx(gd, sample, verbose):
    """returns dict of strand:index"""
    indexes = {}
    # parse the file name
    sample = sample.split(".peaks")[0]
    with Genome(gd) as genome:
        for i, track in enumerate(genome.tracknames_continuous):
            # t_sample, strand, extraneous = track.split(".", 2)
            if "neg" in track:
                t_sample = track.split(".neg")[0]
                strand = "-"
            else:
                t_sample = track.split(".pos")[0]
                strand = "+"
            if sample == t_sample:
                if verbose: sys.stderr.write(">> found track: %s\n" % track)
                indexes[strand] = i
    # should only have one pos and one neg for a sample
    assert(len(indexes.keys()) == 2)
    return indexes

def main(args):
    # pytables via genomedata raises warn PerformanceWarning
    if not args.warnings:
        warnings.simplefilter("ignore")
    track_idx = get_track_idx(args.genomedata, os.path.basename(args.peaks), args.verbose)
    kwargs = {'bed':args.peaks,
              'genomedata':args.genomedata,
              'track_idx':track_idx,
              'peak_ext':args.width / 2,
              'verbose':args.verbose}
    trim_peaks(**kwargs)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('peaks', help='sorted peaks file (bed)')
    p.add_argument('genomedata', help="genome data archive")
    p.add_argument("-w", "--width", default=50, type=int,
                    help="full peak width [ %(default)s ]")
    p.add_argument("-v", "--verbose", action='store_true',
                    help="maximum verbosity [ %(default)s ]")
    p.add_argument("--warnings", action='store_true',
                    help="show warnings [ %(default)s ]")
    main(p.parse_args())