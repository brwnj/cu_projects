#!/usr/bin/env python
# encoding: utf-8
"""
Create pooled count files for testing.
"""
import os
import sys
import gzip
import argparse
from itertools import izip
from toolshed import reader
from collections import defaultdict

def main(count_files, metadata):
    pools = defaultdict(list)
    for toks in reader(metadata):
        for k, v in toks.iteritems():
            if k.startswith("Pool") and v == "TRUE":
                # get the samples
                pool_name = k.split("_")[-1]
                pools[pool_name].append(toks['alias'])
    for pool, samples in pools.iteritems():
        print >>sys.stderr, ">> processing", pool
        for strand in ["pos", "neg"]:
            files = [f for f in count_files if os.path.basename(f).split(".")[0] in samples and strand in os.path.basename(f)]
            pool_file = gzip.open("{pool}.{strand}.txt.gz".format(pool=pool, strand=strand), "wb")
            for count_data in izip(*[reader(f, header=False) for f in files]):
                total_count = sum([int(c[-1]) for c in count_data])
                print >>pool_file, "{name}\t{count}".format(name=count_data[0][0], count=total_count)
            pool_file.close()

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("metadata")
    p.add_argument("counts", nargs="+")
    args = p.parse_args()
    main(args.counts, args.metadata)
