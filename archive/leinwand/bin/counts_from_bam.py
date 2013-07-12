#! /usr/bin/env python
# encoding: utf-8
"""counts of mRNAs per sample where the reads have been aligned to mirbase."""
import sys
import pybedtools
import pandas as pd
import os.path as op
from collections import Counter

def main(args):
    samples = {}
    for bamfile in args.bam:
        # parse file name for sample id
        sample_name = op.basename(bamfile).split(".")[0]
        if args.verbose:
            sys.stderr.write("Processing %s\n" % sample_name)
        samples[sample_name] = Counter()
        for record in pybedtools.BedTool(bamfile):
            # counts per mRNA per sample
            samples[sample_name].update([record.chrom])
    # create dataframe
    df = pd.DataFrame(samples)
    df.to_csv(sys.stdout, sep="\t")

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("bam", metavar="BAM", nargs="+", help="bam or bams")
    p.add_argument("-v", "--verbose", action="store_true",
            help="maximum verbosity")
    main(p.parse_args())