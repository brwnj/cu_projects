#!/usr/bin/env python
# encoding: utf-8
"""
From notSIF format, outputs annotated Cytoscape edge attribute file.
"""

import argparse
import sys
import tempfile
from pybedtools import BedTool
from toolshed import reader
from itertools import izip

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def annotate_peaks(notsif, beds, names):
    """Takes notsif, transforms to bed, and outputs annotation of where the 
    miRNA seed is interrogating via Cytoscape edge attribute file.
    """
    strand = find_strand_from_filename(notsif)
    
    mirna_bed = BedTool(notsif_to_bed(notsif, strand), from_string=True)

    # create the reference beds
    reference = {}
    for name, bed in izip(names, beds):
        reference[name] = BedTool(bed)
    
    for name in names:
        
        # intersect the mirna bed with the reference annotations
        for hit in mirna_bed.intersect(reference[name], s=True, stream=True):
            
            # name field returned from notsif_to_bed is delimited by "|"
            mirna_name = hit.name.split("|")[0]
            gene_name = hit.name.split("|")[1]
            # Cytoscape formatting
            seed_length = '(%s)' % hit.score
            fields = (mirna_name, seed_length, gene_name, "=", name)
            print " ".join(map(str, fields))


def notsif_to_bed(notsif, strand):
    """Converts notsif into bed with miRNA and gene both in the name."""
    # tempfile = 'tmpbed.bed'
    bed = ""
    for tok in reader(notsif, header="name chrom start stop gene".split()):
        # length = seed stop position minus 1
        length = int(tok['name'].rsplit("-", 1)[1]) - 1
        mirna_name = tok['name'].split("|")[0]
        
        bedline = "%s\t%s\t%s\t%s|%s\t%d\t%s\n" % (tok['chrom'], \
                        tok['start'], tok['stop'], mirna_name, tok['gene'], \
                        length, strand)
        # i hope no one ever sees this.
        bed = bed + bedline
    return bed


def find_strand_from_filename(filename):
    if 'pos' in filename:
        return '+'
    elif 'neg' in filename:
        return '-'
    else:
        raise TypeError, "Cannot learn strand from bed filename"


def main():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('notsif', help='cliptools-aggregate-peaks output')
    p.add_argument('--beds', nargs='+', 
        help='reference beds for each feature to annotate')
    p.add_argument('--names', nargs='+', 
        help='names corresponding to reference files')
    args = p.parse_args()
    if not args.names or not args.beds:
        sys.exit(p.print_help())
    annotate_peaks(args.notsif, args.beds, args.names)


if __name__ == "__main__":
    main()