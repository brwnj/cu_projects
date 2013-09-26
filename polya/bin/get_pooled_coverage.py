#!/usr/bin/env python
# encoding: utf-8
"""
get pooled coverage...

i'm doing this after creating the script to get the counts. doing this first
would have enabled getting the counts in a much simpler fashion.
"""

import os
import sys
import argparse
from bsub import bsub
from toolshed import reader
from collections import defaultdict

def main(bams, chromsizes, metadata, projectid):
    pools = defaultdict(list)
    for toks in reader(metadata):
        for k, v in toks.iteritems():
            if k.startswith("Pool") and v == "TRUE":
                # get the samples
                pool_name = k.split("_")[-1]
                pools[pool_name].append(toks['alias'])
    print pools
    merged_bams = []
    jobids = []
    submit = bsub("mergebam", P=projectid)
    for pool, samples in pools.iteritems():
        print >>sys.stderr, ">> processing", pool
        # merge the bams
        bam_files = [bam for bam in bams if os.path.basename(bam).split(".")[0] in samples and "UMIs_not_removed" not in bam and "pos" not in bam and "neg" not in bam]
        merged_bam = "{pool}.bam".format(pool=pool)
        merged_bams.append(merged_bam)
        cmd = "samtools merge {merged_bam} ".format(merged_bam=merged_bam) + " ".join(bam_files)
        print cmd
        # jobids.append(submit(cmd))
        # convert bam using bam2bw
        
        # bam2bw.py -5 -b -v {umibam} {chrom_sizes} {project_id}

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("metadata")
    p.add_argument("chromsizes")
    p.add_argument("bams", nargs="+")
    p.add_argument("--project_id", dest="projectid", default="None")
    args = p.parse_args()
    main(args.bams, args.chromsizes, args.metadata, args.projectid)
