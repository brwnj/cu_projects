#!/usr/bin/env python
# encoding: utf-8
"""
Takes miRNA coords in BED format. Outputs:

hsa-miR-429 = 11872.0
hsa-miR-4251 = 0
hsa-miR-4689 = 5.0
"""

import argparse
import os
import sys
from toolshed import reader
from pybedtools import BedTool
from genomedata import Genome
from numpy import nansum
from collections import defaultdict


def peak_intensity(bed, genomedata, track, sif):
    """For each BED entry (miRNA coordinate), get the max peak intensity.
    
    bed - reference file with miRNA name, coordinate, and strand
    genomedata - intensity archive
    track - information organized in the archive as tracks
    sif - network file. anything outside of this should be ignored
    """
    bed = BedTool(bed)
    
    network = mirs_in_sif(sif)
    mir_tracker = defaultdict(list)
    with Genome(genomedata) as genome:
        for b in bed:
            try:
                # filter out any chromosomes not in the genome data archive
                chromosome = genome[b.chrom]
                # filter out any mirs not present in the network
                lookup = network[b.name]
            except KeyError:
                continue
            
            strand = ""
            if b.strand is "+":
                strand = "pos"
            else:
                strand = "neg"
            
            strand_corrected = []
            for i, t in enumerate(track):
                strand_corrected.append("%s.%s" % (t, strand))
            
            # intensities across specified tracks for peak region
            intensities = chromosome[b.start:b.stop, strand_corrected]
            max_intensity = get_peak_max(intensities)
            
            mir_tracker[b.name].append(max_intensity)
            
            # this was for printing bed format
            #fields = (b.chrom, b.start, b.stop, b.name, max_intensity, b.strand)
            #print "\t".join(map(str, fields))
    
    # taking only the maximum if a mir has multiple locations
    for mirname, intensities in mir_tracker.iteritems():
        fields = (mirname, "=", float(max(intensities)))
        print " ".join(map(str, fields))


def get_peak_max(intensities):
    """Finds the maximum peak intensity value from ndarray.
    """
    maximum = 0
    for intensity in intensities:
        # sum intensity across samples
        sum_intensity = nansum(intensity)
        
        if sum_intensity > maximum:
            maximum = sum_intensity
    return maximum


def mirs_in_sif(sif):
    """returns dictionary of miRNAs present in the network."""
    network = {}
    for s in reader(sif, header="mir op gene".split()):
        network[s['mir']] = s['gene']
    return network


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    
    p.add_argument('bed', help='miRNA coordinates')
    p.add_argument('sif', help='Cytoscape network file')
    p.add_argument('genomedata', help="genome data archive")
    p.add_argument('trackname', nargs='+',
                    help='the names of tracks corresponding to study case')

    args = p.parse_args()
    peak_intensity(args.bed, args.genomedata, args.trackname, args.sif)


if __name__ == "__main__":
    main()