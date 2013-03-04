#!/usr/bin/env python
# encoding: utf-8
"""
find regions of enrichment. doesn't care about joining adjacent windows since
that will be completed later using comb-p by brent p.

todo: minimum summit height to be considered
"""
import sys
import warnings
from scipy.stats import poisson
from genomedata import Genome
from numpy import nansum, isnan

def call_peaks(genomedata, track, width, p_cutoff, verbose):
    with Genome(genomedata) as genome:
        strand = find_strand(track)
        track_idx = genome.index_continuous(track)
        # effective genome size
        size = sum([chrom.end for chrom in genome])
        # coverage across all chromosome for one track
        coverage = genome.sums[track_idx] / size
        # global lambda
        g_lambda = coverage * width
        if verbose:
            sys.stderr.write(">> effective genome size: %f\n" % size)
            sys.stderr.write(">> observed coverage: %f\n" % coverage)
            sys.stderr.write(">> calculated global lambda: %f\n" % g_lambda)
        # need this for unique peak names throughout genome
        i = 0
        for chromosome in genome:
            if verbose: sys.stderr.write(">> processing %s\n" % chromosome.name)
            for supercontig, continuous in chromosome.itercontinuous():
                whole_track = continuous[:, track_idx]
                # possible peaks through `whole_track`
                starts = xrange(0, len(whole_track), width)
                for start in starts:
                    stop = start + width
                    # location is relative within this continous
                    counts = whole_track[start:stop]
                    # skip where no reads
                    if isnan(nansum(counts)): continue
                    # calculate p-value
                    p = poisson.pmf(nansum(counts), g_lambda)
                    # ignore garbage
                    if p > p_cutoff: continue
                    i += 1
                    # print the current peak, fixing location
                    fields = [chromosome.name, supercontig.start + start,
                                supercontig.start + stop, "%s_%d" % (track, i),
                                p, strand]
                    print "\t".join(map(str, fields))

def find_strand(track):
    """determine strand symbol from track name."""
    strand = "."
    if "pos" in track: strand = "+"
    if "neg" in track: strand = "-"
    return strand        

def main(args):
    kwargs = {'genomedata':args.genomedata,
              'track':args.track,
              'width':args.width,
              'p_cutoff':args.p_value,
              'verbose':args.verbose}
    # pytables via genomedata warns PerformanceWarning
    if not args.warnings:
        warnings.simplefilter("ignore")
    call_peaks(**kwargs)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('genomedata', help="genomedata archive")
    p.add_argument('track', help='name of the track to call peaks across')
    p.add_argument('-p', '--p-value', dest="p_value", default=0.00001, type=int, 
                    help="ignore peaks above this cutoff [ %(default)s ]")
    p.add_argument('-w', '--width', default=50, type=int,
                    help="window size or close to expected peak width \
                    [ %(default)s ]")
    p.add_argument("-v", "--verbose", action='store_true',
                    help="maximum verbosity [ %(default)s ]")
    p.add_argument("--warnings", action='store_true',
                    help="show warnings [ %(default)s ]")
    main(p.parse_args())