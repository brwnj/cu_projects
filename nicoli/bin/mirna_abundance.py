#!/usr/bin/env python
# encoding: utf-8
"""
Takes miRNA coords in BED format. Outputs:

chrom, start, stop, miRNA, score, strand, miRNA_max_abundance
"""
import sys
from numpy import nansum
from toolshed import reader
from genomedata import Genome
from collections import defaultdict
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def peak_max(intensities):
    """Finds the maximum peak intensity value from ndarray.
    """
    maximum = 0
    for intensity in intensities:
        # sum intensity across samples
        sum_intensity = nansum(intensity)

        if sum_intensity > maximum:
            maximum = sum_intensity
    return maximum

def present_mirnas(files):
    """returns set of observed miRNAs given list of seed mapping files.
    """
    seen = set()
    for f in files:
        for toks in reader(f, header=['mirna', 'chrom', 'start', 'stop', 'gene']):
            mirna_name = toks['mirna'].split("|", 1)[0]
            seen.add(mirna_name)
    return seen

def main(bed, seed_mapping_pos, seed_mapping_neg, genomedata, trackname_pos, trackname_neg):
    """For each BED entry (miRNA coordinate), get the max peak intensity.
    """
    bed_header = ['chrom','start','stop','name','score','strand']
    mirna_reference = dict()
    observed_mirnas = present_mirnas([seed_mapping_pos, seed_mapping_neg])
    mirna_abundances = defaultdict(list)
    with Genome(genomedata) as genome:
        for chrom in genome:
            for b in reader(bed, header=bed_header):
                # faster to load chromosome of genomedata archive
                if b['chrom'] != chrom.name: continue
                # strand info needed for genomedata lookup
                assert b['strand'] in {"+", "-"}
                # only interested in seeds with hits
                if not b['name'] in observed_mirnas: continue

                start = int(b['start'])
                stop = int(b['stop'])
                mirna_reference[b['name']] = "\t".join(b[tok] for tok in bed_header)

                if b['strand'] == "+":
                    intensities = chrom[start:stop, trackname_pos]
                else:
                    intensities = chrom[start:stop, trackname_neg]
                max_intensity = peak_max(intensities)

                mirna_abundances[b['name']].append(max_intensity)

    # print bed6+
    for mirna, intensities in mirna_abundances.iteritems():
        max_abundance = float(max(intensities))
        if max_abundance == 0: continue
        print "{bed6}\t{max_abundance}".format(bed6=mirna_reference[mirna],
                    max_abundance=max_abundance)

if __name__ == "__main__":
    p = ArgumentParser(description=__doc__,
                   formatter_class=RawDescriptionHelpFormatter)
    p.add_argument('bed', help='miRNA reference')
    p.add_argument('seed_mapping_pos', help='output of `cliptools-aggregate-peaks`')
    p.add_argument('seed_mapping_neg', help='output of `cliptools-aggregate-peaks`')
    p.add_argument('genomedata', help="genome data archive")
    p.add_argument('trackname_pos', help='corresponding positive strand track name in genomedata archive')
    p.add_argument('trackname_neg', help='corresponding positive strand track name in genomedata archive')

    args = vars(p.parse_args())
    main(**args)
