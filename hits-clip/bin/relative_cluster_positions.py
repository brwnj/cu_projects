#!/usr/bin/env python
# encoding: utf-8
"""
Prints distance, as a fraction of peak position within any given transcript,
of Ago-mRNA peaks relative to a feature
"""
import argparse
import sys
from pybedtools import BedTool
from collections import defaultdict

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


class ComplexName(object):
    """Expecting UCSC mRNA bed with name field like:
    
    AK127011_exon_0_0_chr1_67051162_r
    
    """
    def __init__(self, longname):
        super(ComplexName, self).__init__()
        attrs = longname.rsplit("_",6)
        self.name = attrs[0]
        self.feature = attrs[1]
        self.number = int(attrs[2])


class ComplexLine(object):
    """After all the bedtools actions are concatenated onto the original."""
    def __init__(self, interval):
        super(ComplexLine, self).__init__()
        self.peakchrom = interval[0]
        self.peakstart = int(interval[1])
        self.peakstop = int(interval[2])
        self.peakname = interval[3]
        self.peakstrand = interval[4]
        self.exonchrom = interval[5]
        self.exonstart = int(interval[6])
        self.exonstop = int(interval[7])
        self.exoninfo = ComplexName(interval[8])
        self.exonscore = interval[9]
        self.exonstrand = interval[10]
        self.exonoverlap = interval[11]
        self.gtfchrom = interval[12]
        self.gtfsource = interval[13]
        self.gtffeature = interval[14]
        self.gtfstart = int(interval[15])
        self.gtfstop = int(interval[16])
        self.gtfscore = interval[17]
        self.gtfstrand = interval[18]
        self.gtfframe = interval[19]
        self.gtfattrs = interval[20]
        self.gtfdistance = int(interval[21])
        # nothing returned from `closest` for this peak
        if self.gtfattrs == ".":
            raise ValueError


def make_exon_lib(bed):
    """Given bed, converts to dictionary of refseq name and exon lengths.
    
    returns defaultdict(dict)
    """
    bed = BedTool(bed)
    exon_lib = defaultdict(dict)
    for b in bed:
        try:
            exon = ComplexName(b.name)
            exon_lib[exon.name][exon.number] = b.length
        except ValueError:
            # omit chr??_random
            pass
    return exon_lib


def main():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('peaks', help='peaks bed')
    p.add_argument('exons', help='refseq exons from UCSC')
    p.add_argument('gtf', help='refseq gtf with feature of interest')
    p.add_argument('feature', help='feature of interest in the gtf')
    p.add_argument('-v', '--verbose', action="store_true", help='maximum verbosity')
    args = p.parse_args()
    
    if args.verbose: sys.stderr.write(">> building exon library...\n")
    exon_lib = make_exon_lib(args.exons)
    
    peaks = BedTool(args.peaks)
    exons = BedTool(args.exons)
    full_ref = BedTool(args.gtf)
    
    if args.verbose: sys.stderr.write(">> filtering for feature...\n")
    filtered_ref = full_ref.filter(lambda gtf: gtf[2] == args.feature)
    
    if args.verbose: sys.stderr.write(">> selecting exonic peaks...\n")
    exonic_peaks = peaks.intersect(exons, wo=True)
    
    if args.verbose: sys.stderr.write(">> calculating distance fractions...\n")
    # D for distance (returns negative if upstream)
    for peak in exonic_peaks.closest(filtered_ref, D="a"):
        try:
            p = ComplexLine(peak)
            corrected_distance = 0.0
            total_exon_length = 0.0
            # parse gtf attrs
            gene_id = p.gtfattrs.split(';')[0].rstrip('"').lstrip('gene_id "')

            # looking downstream wrt peak
            if p.gtfdistance > 0:
                # exon with peak
                corrected_distance = p.exonstop - p.peakstop
                for exon in exon_lib[p.exoninfo.name]:
                    # add downstream exon lengths
                    if exon > p.exoninfo.number:
                        corrected_distance += exon_lib[p.exoninfo.name][exon]
                        
            # looking upstream wrt peak
            else:
                # exon with peak
                corrected_distance = p.peakstart - p.exonstart
                for exon in exon_lib[p.exoninfo.name]:
                    # add upstream exon lengths
                    if exon < p.exoninfo.number:
                        corrected_distance += exon_lib[p.exoninfo.name][exon]
            
            for exon in exon_lib[p.exoninfo.name]:
                total_exon_length += exon_lib[p.exoninfo.name][exon]
            
            # fraction
            print (corrected_distance / total_exon_length)
        
        except ValueError:
            continue


if __name__ == "__main__":
    main()