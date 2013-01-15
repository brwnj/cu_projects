#!/usr/bin/env python
# encoding: utf-8
"""
parse novoalign alignment output to produce statistics table.
"""
import sys
import pandas as pd
import os.path as op
from toolshed import nopen

def parse_stats(fin):
    """parse novoalign stderr alignment metrics. returns dict of metrics."""
    stats = {}
    sample = ""
    # get the sample name
    for line in nopen(fin):
        line = line.rstrip("\n\r").split()
        if 'novoalign' in line[1]:
            for i, flag in enumerate(line):
                if "-f" in flag:
                    # assumes sample is first in file name, delimited by periods
                    sample = op.basename(line[i + 1]).split(".")[0]
                    break
    # get the data
    for line in nopen(fin):
        line = line.rstrip("\n\r").split()
        if line[1] == "Read":
            stats["total_reads"] = line[-1]
        if line[1] == "Aligned:":
            stats["aligned"] = line[-1]
        if line[1] == "Unique":
            stats["uniquely_aligned"] = line[-1]
    return sample, stats

def main(args):
    stats = {}
    for f in args.files:
        sample_name, mapping_stats = parse_stats(f)
        stats[sample_name] = mapping_stats
        mapping_stats['percent_aligned'] = "%0.2f" % \
            (float(mapping_stats['aligned']) / \
                int(mapping_stats['total_reads']) * 100)
        mapping_stats['percent_unique'] = "%0.2f" % \
            (float(mapping_stats['uniquely_aligned']) / \
                int(mapping_stats['total_reads']) * 100)
    df = pd.DataFrame(stats)
    df = df.ix[["total_reads","aligned","percent_aligned","uniquely_aligned","percent_unique"]]
    df.to_csv(sys.stdout, sep="\t")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", metavar="FILES", nargs="+", 
            help="stderr of novoalign mapping")
    args = p.parse_args()
    main(args)