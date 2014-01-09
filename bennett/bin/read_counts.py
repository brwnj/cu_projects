#!/usr/bin/env python
# encoding: utf-8
"""
counts of reads per barcode
counts of reads after umi filtering
counts of reads after trimming adapter sequences
counts of reads after joining
"""
import sys
import pandas
import argparse
from toolshed import nopen

def flen(fname):
    return sum(1 for line in nopen(fname)) / 4

def gsample(fname):
    return fname.split(".", 1)[0]

def progress(position, total):
    print >>sys.stderr, "  >> processing %d of %d" % (position, total)

def main(unprocessed_fastqs, umi_filtered_fastqs, adapter_trimmed_fastqs, joined_fastqs):
    assert len(unprocessed_fastqs) == len(umi_filtered_fastqs) == len(adapter_trimmed_fastqs) == len(joined_fastqs)
    to_process = len(unprocessed_fastqs)
    out_order = ["read count", "umi filtering", "adapter trimmed", "joined"]
    metadata = dict()
    print >>sys.stderr, ">> reading unprocessed fastqs"
    for i, fq in enumerate(unprocessed_fastqs, start=1):
        progress(i, to_process)
        sample = gsample(fq)
        metadata[sample] = {"read count":flen(fq)}
    print >>sys.stderr, ">> reading umi filtered fastqs"
    for i, fq in enumerate(umi_filtered_fastqs, start=1):
        progress(i, to_process)
        sample = gsample(fq)
        metadata[sample]["umi filtering"] = flen(fq)
    print >>sys.stderr, ">> reading adapter trimmed fastqs"
    for i, fq in enumerate(adapter_trimmed_fastqs, start=1):
        progress(i, to_process)
        sample = gsample(fq)
        metadata[sample]["adapter trimmed"] = flen(fq)
    print >>sys.stderr, ">> reading joined fastqs"
    for i, fq in enumerate(joined_fastqs, start=1):
        progress(i, to_process)
        sample = gsample(fq)
        metadata[sample]["joined"] = flen(fq)
    df = pandas.DataFrame(metadata)
    df.loc[out_order,:].to_csv(sys.stdout, sep="\t")

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('--one', nargs="+", help="unprocessed fastqs")
    p.add_argument('--two', nargs="+", help="umi filtered fastqs")
    p.add_argument('--three', nargs="+", help="biological adapter trimmed")
    p.add_argument('--four', nargs="+", help="joined reads fastqs")
    args = p.parse_args()
    main(args.one, args.two, args.three, args.four)
