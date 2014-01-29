#!/usr/bin/env python
# encoding: utf-8
"""
chr1	11868	14412	TCONS_00000004	0	+	11868	14412	0	5	189,49,109,150,960,	0,310,744,1356,1584,
chr1	11868	14412	TCONS_00000003	0	+	11868	14412	0	5	189,49,127,253,752,	0,310,726,1534,1792,
chr1	11868	14409	TCONS_00000001	0	+	11868	14409	0	3	359,109,1189,	0,744,1352,
"""
import re
import sys
import argparse
from toolshed import reader
from itertools import groupby

class Bed12(object):
    """Convert BedTool interval to Bed12."""
    __slots__ = ["chrom", "start", "end", "name", "score", "strand",
                    "thickStart", "thickEnd", "itemRgb", "blockCount",
                    "blockSizes", "blockStarts"]
    def __init__(self, gtf, attr):
        self.chrom = gtf.chrom
        self.start = gtf.start - 1
        self.end = gtf.end
        self.name = attr
        self.score = gtf.score
        self.strand = gtf.strand
        self.thickStart = self.start
        self.thickEnd = self.end
        self.itemRgb = 0
        self.blockCount = 1
        self.blockSizes = "{size},".format(size=gtf.end - gtf.start)
        self.blockStarts = "{start},".format(start=0)

    def __repr__(self):
        return "Bed12({chrom}:{start}-{end})".format(chrom=self.chrom,
                start=self.start, end=self.end)
    
    def __str__(self):
        return "\t".join(str(getattr(self, s)) for s in self.__slots__)
    
    def add_exon(self, gtf, attr):
        self.end = gtf.end
        self.thickEnd = gtf.end
        self.blockSizes += "{size},".format(size=gtf.end - gtf.start)
        self.blockStarts += "{start},".format(start=gtf.start - self.start)
        self.blockCount += 1

class GTF(object):
    """
    http://uswest.ensembl.org/info/website/upload/gff.html
    """
    __slots__ = ['seqname', 'source', 'feature', 'start', 'end', 'score',
                    'strand', 'frame', 'attributes', 'chrom', 'stop']

    def __init__(self, args):
        for s, v in zip(self.__slots__[:9], args):
            setattr(self, s, v)
        self.start = int(self.start)
        self.end = int(self.end)
        self.chrom = self.seqname
        self.stop = self.end
    
    def __repr__(self):
        return "GTF({seqname}:{start}-{end})".format(seqname=self.seqname,
                start=self.start, end=self.end)

def main(gtf, attribute_id):
    parse_attr = re.compile('{attribute_id} "([^"]+)"'.format(attribute_id=attribute_id))
    for k, group in groupby(reader(gtf, header=GTF), lambda x: parse_attr.findall(x.attributes)[0]):
        # loop over the group building the single bed12 line
        bed = None
        for gtf_entry in group:
            if not bed:
                bed = Bed12(gtf_entry, k)
            else:
                bed.add_exon(gtf_entry, k)
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
