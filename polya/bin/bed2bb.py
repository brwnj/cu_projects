#!/usr/bin/env python
# encoding: utf-8
"""
Convert bed to bigbed. Score field will be scaled to 1000 as max. UCSC
bedToBigBed must be in PATH.
"""
import os
import sys
import argparse
import tempfile
import subprocess as sp
from toolshed import reader

class Bed6(object):
    __slots__ = ['chrom','start','stop','name', 'score', 'strand']
    
    def __init__(self, args):
        for k, v in zip(self.__slots__, args):
            setattr(self, k, v)
        self.score = int(self.score)
    
    def __str__(self):    
        return "\t".join([getattr(self, s) for s in self.__slots__])

def scale(score, omax, omin, smax, smin):
    """
    >>> scale(2871, 4871, 0, 1000, 0)
    589
    """
    return ((smax - smin) * (score - omin) / (omax - omin)) + smin

def main(beds, sizes, mx, mn):
    for bed in beds:
        name = bed.split(".bed")[0]
        scores = []
        # find max and min
        for b in reader(bed, header=Bed6):
            scores.append(b.score)
        try:
            scoremx = float(max(scores))
        except ValueError:
            print >>sys.stderr, bed, "appears empty."
            sys.exit(1)
        scoremn = float(min(scores))
        tmp = open(tempfile.mkstemp(suffix=".bed")[1], 'wb')
        # scale the values
        for b in reader(bed, header=Bed6):
            b.score = str(int(scale(b.score, scoremx, scoremn, mx, mn)))
            print >>tmp, b
        tmp.close()
        #convert to bb
        cmd = "bedToBigBed -type=bed6 {tmp} {sizes} {name}.bb".format(
                tmp=tmp.name, sizes=sizes, name=name)
        sp.call(cmd, shell=True)
        os.remove(tmp.name)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("bed", metavar="BED", nargs="+",
            help="bed6 to convert to bigbed")
    p.add_argument("sizes", metavar="SIZES", help="genome sizes file")
    p.add_argument("--score-max", dest="scoremx", default=1000, type=int,
            help="maximum score to be observed")
    p.add_argument("--score-min", dest="scoremn", default=0, type=int,
            help="mininum score to be observed")
    args = p.parse_args()
    if 0 < args.scoremx > 1000 or 0 > args.scoremn > 1000 \
            or args.scoremn > args.scoremx:
        print >>sys.stderr, "Scores need to satisfy: 0 <= input_val <= 1000"
        sys.exit(1)
    main(args.bed, args.sizes, args.scoremx, args.scoremn)
