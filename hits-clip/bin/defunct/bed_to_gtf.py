#!/usr/bin/env python
# encoding: utf-8
"""
Convert beds to gff format for annotation purposes.

Author: Joe Brown
Date: 2012-05-11

"""
import sys
from pybedtools import BedTool


def features_from_filename(filename):
    """returns the feature name assuming the file name is delimited by periods
    and the feature name is in position 1, e.g. refseq.cds.bed.gz
    """
    parts = filename.split(".")
    source = parts[0]
    feature = parts[1]
    if len(parts) < 3:
        raise TypeError, "File name must not be delimited by periods."
    return source, feature


def main(beds):
    for bed in beds:
        source, feature = features_from_filename(bed)
        bed = BedTool(bed)
        
        if feature.lower() == 'cds':
            frame = "0"
        else:
            frame = "."
        
        for b in bed:
            fields = (b.chrom, source, feature, b.start, b.stop, \
                        b.score, b.strand, frame, b.name)
            print "\t".join(map(str, fields))


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter,
            )
    p.add_argument('bed', nargs='+', help='bed or beds to use.')
    args = p.parse_args()
    
    main(args.bed)