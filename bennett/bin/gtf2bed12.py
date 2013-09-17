#!/usr/bin/env python
# encoding: utf-8
"""
chr1	11868	14412	TCONS_00000004	0	+	11868	14412	0	5	189,49,109,150,960,	0,310,744,1356,1584,
chr1	11868	14412	TCONS_00000003	0	+	11868	14412	0	5	189,49,127,253,752,	0,310,726,1534,1792,
chr1	11868	14409	TCONS_00000001	0	+	11868	14409	0	3	359,109,1189,	0,744,1352,
"""
import sys
import argparse
from itertools import groupby
from pybedtools import BedTool

class Bed12(object):
    """Convert BedTool interval to Bed12."""
    __slots__ = ["chrom", "start", "end", "name", "score", "strand",
                    "thickStart", "thickEnd", "itemRgb", "blockCount",
                    "blockSizes", "blockStarts"]
    def __init__(self, bedtool, attr):
        self.chrom = bedtool.chrom
        self.start = bedtool.start
        self.end = bedtool.end
        self.name = bedtool.attrs[attr]
        self.score = bedtool.score
        self.strand = bedtool.strand
        self.thickStart = bedtool.start
        self.thickEnd = bedtool.end
        self.itemRgb = 0
        self.blockCount = 1
        self.blockSizes = "{size},".format(size=bedtool.end - bedtool.start)
        self.blockStarts = "{start},".format(start=0)

    def __repr__(self):
        return "Bed12({chrom}:{start}-{end})".format(chrom=self.chrom,
                start=self.start, end=self.end)
    
    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)
    
    def add_exon(self, bedtool, attr):
        self.end = bedtool.end
        self.thickEnd = bedtool.end
        self.blockSizes += "{size},".format(size=bedtool.end - bedtool.start)
        self.blockStarts += "{start},".format(start=bedtool.start - self.start)
        self.blockCount += 1

def parse_attrs(attrs):
    attrs = attrs.replace('"', '')
    spaced_attrs = [a.strip() for a in attrs.split(";") if a]
    attributes = {}
    attributes.up
    for k in spaced_attrs:
        key, value = k.split(" ")
        attributes[key] = value
    return attributes
    
    
    [dict(zip(entry.split(" "))) for entry in spaced_attrs]

    [ (v-1, v+1) for v in [ 2*i for i in list ] if someCondition(v) ]
    
    [dict(zip([v[0]],[v[1]])) for v in [entry.split(" ") for entry in spaced_attrs]]

def main(gtf, attr):
    gtf = BedTool(gtf)
    for k, group in groupby(gtf, lambda x: x.attrs[attr]):
        # loop over the group building the single bed12 line
        gtf_entries = sorted([rec for rec in group])
        
        for gtf_entry in gtf_entries:
            
            
            if not bed:
                bed = Bed12(gtf_entry, attr)
            else:
                bed.add_exon(gtf_entry, attr)
        print bed

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('gtf', help="gtf file sorted and grouped by `attr`. if it's \
            the result of a cufflinks tool, it should already be sorted in the \
            proper order.")
    p.add_argument('attr', help="attribute ID by which to group bed12 entries.")
    p.add_argument('--feature-type', default='exon', help="feature type to \
            join -- all others are filtered out")
    args = p.parse_args()
    main(args.gtf, args.attr)
