#!/usr/bin/env python
# encoding: utf-8
"""
Writes an edge attributes file with peak intensity for each miRNA seed match.
Requires the file output from ``cliptools-aggregrate-peaks``, which I refer
to as 'notsif'.

Input file:
hsa-miR-483-3p|seed:2-8	chrY	22897163	22897223	RBMY1F
mmu-miR-148a*|seed:2-8	chrY	22897163	22897223	RBMY1F

Output file:
PeakIntensity
hsa-miR-483-3p (8) RBMY1F = 35.0
mmu-miR-148a* (8) RBMY1F = 12.0
"""

import argparse
import sys
from toolshed import reader
from genomedata import Genome
from numpy import nansum

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def peak_attributes(notsif, genomedata, track):
    """For each mirna peak match, find the max intensity and write
    Cytoscape attrs file.
    """
    with Genome(genomedata) as genome:
        
        # saving lookups for all the repetition
        previous_chrom = previous_start = \
            previous_stop = previous_intensity = None
        
        for line in reader(notsif, header="name chrom start stop gene".split()):
            # length = seed stop position minus 1
            length = int(line['name'].rsplit("-", 1)[1]) - 1
            
            # format for output
            seed_length = "(%d)" % length
            mirna_name = line['name'].split("|")[0]
            
            if previous_chrom and \
                    line['chrom'] is previous_chrom and \
                    line['start'] is previous_start and \
                    line['stop'] is previous_stop:
                
                # new micro, but interrogating same peak as before
                max_intensity = previous_intensity

                fields = (mirna_name, seed_length, \
                            line['gene'], "=", max_intensity)
                print " ".join(map(str, fields))
                continue
            
            chromosome = genome[line['chrom']]
            
            # intensities across specified tracks for peak region
            intensities = chromosome[int(line['start']):int(line['stop']), \
                                        track]
            max_intensity = get_max(line['chrom'], line['start'], intensities)

            previous_intensity = max_intensity
            previous_chrom = line['chrom']
            previous_start = line['start']
            previous_stop = line['stop']
            
            fields = (mirna_name, seed_length, line['gene'], \
                        "=", max_intensity)
            # output peak attributes file
            print " ".join(map(str, fields))


def get_max(chrom, start, intensities):
    """Finds the maximum intensity value."""
    maximum = 0
    for intensity in intensities:
        
        # sum intensity across samples
        sum_intensity = nansum(intensity)
        if sum_intensity > maximum:
            maximum = sum_intensity
    
    return maximum


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('notsif', help='output of cliptools-aggregate-peaks')
    
    p.add_argument('genomedata', help='genome data archive')
    
    p.add_argument('trackname', nargs='+', 
                    help='track names corresponding to study case')

    args = p.parse_args()
    peak_attributes(args.notsif, args.genomedata, args.trackname)


if __name__ == "__main__":
    main()