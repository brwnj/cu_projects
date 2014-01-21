#!/usr/bin/env python
# encoding: utf-8
"""
Trim peaks by finding the point of maximum intensity for a given interval within
the genomedata archive then expanding to `width`.
"""

import os
import sys
import warnings
from numpy import nansum
from genomedata import Genome
from pybedtools import BedTool
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


def get_mid(start, counts):
    """Finds median position of the maximum intensity value."""

    maximum = 0
    maxpositions = []

    # get all of the positions with the max intensity
    for idx, count in enumerate(counts, start=start):
        summed = nansum(count)
        if summed < maximum: continue
        if summed == maximum:
            maxpositions.append(idx)
            continue
        if summed > maximum:
            maximum = summed
            del maxpositions[0:len(maxpositions)]
            maxpositions.append(idx)

    # return midpoint among the max intensities
    midpoint = int(maxpositions[len(maxpositions) / 2])

    return midpoint


def indexes_from_archive(genomedata, track_names, verbose):
    """returns list of indexes."""
    if verbose:
        print >>sys.stderr, ">> retrieving track indexes"
    indexes = []
    with Genome(gd) as genome:
        for track_name in track_names:
            try:
                idx = genome.index_continuous(track_name)
                indexes.append(idx)
            except KeyError:
                print >>sys.stderr, "invalid track name:", track_name
                sys.exit(1)
    return indexes


def main(peaks, genomedata, tracks, width, verbose, warnings):
    # pytables via genomedata raises warn PerformanceWarning
    if not warnings:
        warnings.simplefilter("ignore")

    track_indexes = indexes_from_archive(genomedata, tracks, verbose)

    # half of width to extend start and stop
    from_start = width / 2
    from_stop = from_start + 1

    bed = BedTool(bed)

    with Genome(genomedata) as gd:
        for chrom in gd:
            if verbose:
                print >>sys.stderr, ">> processing {chrom}...".format(chrom=chrom)

            for supercontig, continuous in chrom.itercontinuous():

                track_continuous = continuous[:, [track_indexes]]

                # will be iterating over this for each supercontig
                for b in bed:
                    # there must be a better way...
                    if b.chrom != chrom.name: continue

                    # only trim within this supercontig
                    if not b.start >= supercontig.start: continue
                    if not b.stop <= supercontig.end: continue

                    # get the relative position since continous will start at 0
                    relative_start = supercontig.project(b.start)
                    relative_stop = supercontig.project(b.stop)

                    # count slice for bed region
                    counts = track_continuous[relative_start:relative_stop]

                    # pass the actual start since we're trimming the peak here
                    mid = get_mid(b.start, counts)
                    start = mid - from_start
                    stop = mid + from_stop
                    print "{chrom}\t{start}\t{stop}\t{name}\t{score}\t{strand}".format(
                            chrom=b.chrom, start=start, stop=stop, name=b.name,
                            score=b.score, strand=b.strand)


if __name__ == "__main__":
    p = ArgumentParser(description=__doc__,
                    formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('peaks', help='sorted peaks as bed')
    p.add_argument('genomedata', help="genomedata archive")
    p.add_argument("-t", "--trackname", action="append", required=True,
                    help=("track names for samples comprising the replicates"
                            "in the peaks bed"))
    p.add_argument("-w", "--width", default=50, type=int,
                    help="full peak width")
    p.add_argument("-v", "--verbose", action='store_true',
                    help="maximum verbosity")
    p.add_argument("--warnings", action='store_true',
                    help="show warnings")
    args = vars(p.parse_args())
    main(**args)
