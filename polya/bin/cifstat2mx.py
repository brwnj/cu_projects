#!/usr/bin/env python
# encoding: utf-8
"""
matrix values are the difference between total zeros and number of clusters
"""
import sys
import argparse
import pandas as pd
from toolshed import nopen
from itertools import islice

class CifStat(object):
    def __init__(self, args):
        self.filepath = args[0].split("/")
        self.ncluster = int(args[1].split()[-1])
        self.zeros = int(args[5].split()[-1])
    
    @property
    def cycle(self):
        return int(self.filepath[3].split(".")[0].lstrip("C"))
    
    @property
    def tile(self):
        return int(self.filepath[-1].split(".")[0].split("_")[-1])

def readstat(cifstat):
    with nopen(cifstat) as fh:
        clean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(clean, 6)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 6
            yield CifStat(rd)

def main(cifstat):
    d = {}
    for r in readstat(cifstat):
        try:
            d[r.cycle][r.tile] = r.ncluster - r.zeros
        except KeyError:
            d[r.cycle] = {}
            d[r.cycle][r.tile] = r.ncluster - r.zeros
    df = pd.DataFrame(d)
    df.to_csv(sys.stdout, sep="\t", index_label="tile")

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('cifstat', help="cifstat with echoed info")
    args = p.parse_args()
    main(args.cifstat)
