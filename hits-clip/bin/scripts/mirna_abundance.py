#!/usr/bin/env python
# encoding: utf-8
"""
Calculate miRNA abundances for singletons or groups with replicates.
"""


import os
import sys
import argparse
from numpy import nansum
from toolshed import reader
from genomedata import Genome
from pybedtools import BedTool
from collections import defaultdict


def peak_intensity(bed, genomedata, trackname, txt):
    """For each BED entry (miRNA coordinate), get the max peak intensity.

    bed - reference file with miRNA name, coordinate, and strand
    genomedata - genomedata archive
    trackname - names of tracks for this sample or group
    txt - mapped seeds and their gene location
    """
    mir_bed = BedTool(bed)

    observed_mirs = parse_txt(txt)
    observed_intensities = defaultdict(list)
    with Genome(genomedata) as genome:
        for b in mir_bed:
            # filter out any mirs not present in the network
            if not b.name in observed_mirs: continue

            chromosome = genome[b.chrom]

            strand = "pos" if b.strand == "+" else "neg"
            selected_tracks = [t for t in trackname if strand in t]

            # intensities across specified tracks for peak region
            intensities = chromosome[b.start:b.stop, selected_tracks]
            max_intensity = get_peak_max(intensities)
            observed_intensities[b.name].append(max_intensity)

    # taking only the maximum if a mir has multiple locations
    for mir_name, values in observed_intensities.iteritems():
        # values is a list of what was observed for a given miRNA
        # by name. some are listed in the bed more than once,
        # so just take the max from whichever region
        print "{name}\t{max_value}".format(name=mir_name, max_value=max(values))


def get_peak_max(intensities):
    """Finds the maximum peak intensity value from ndarray.
    """
    sums = [nansum(i) for i in intensities]
    return max(sums)


def parse_txt(txt):
    """returns dictionary of miRNAs present in the network."""
    observed_mirs = set()
    for t in txt:
        for toks in reader(t, header=['name', 'chrom', 'start', 'stop', 'gene']):
            mir_name = toks['name'].split("|")[0]
            observed_mirs.add(mir_name)
    return observed_mirs


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument('bed', help='miRNA coordinates as bed file')
    p.add_argument('genomedata', help="genome data archive")
    p.add_argument('--tracks', nargs='+', help='the names of tracks corresponding to sample; including positive and negative tracks')
    p.add_argument('--txt', nargs="+", help='mapped seeds text file(s). format: mir|seed, chr, start, stop, gene.')

    args = p.parse_args()
    peak_intensity(args.bed, args.genomedata, args.tracks, args.txt)


if __name__ == "__main__":
    main()
