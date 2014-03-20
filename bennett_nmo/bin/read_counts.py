#!/usr/bin/env python
# encoding: utf-8
"""
counts of reads per barcode
counts of reads after umi filtering
counts of reads after trimming adapter sequences
counts of reads after joining
"""
import os
import sys
import pandas
import argparse
from toolshed import nopen
from collections import defaultdict

def flen(fname):
    return sum(1 for line in nopen(fname)) / 4

def gsample(fname):
    return os.path.basename(fname).split(".", 1)[0].split("_R1", 1)[0]

def progress(position, total):
    print >>sys.stderr, "  >> processing %d of %d" % (position, total)

def process_samples(files, file_type, metadata, key):
    print >>sys.stderr, ">> processing", file_type
    to_process = len(files)
    for i, fq in enumerate(files, start=1):
        progress(i, to_process)
        sample = gsample(fq)
        metadata[sample][key] = flen(fq)

def main(unprocessed_fastqs, umi_filtered_fastqs, adapter_trimmed_fastqs, joined_fastqs):
    assert len(unprocessed_fastqs) == len(umi_filtered_fastqs) == len(adapter_trimmed_fastqs) == len(joined_fastqs)
    out_order = ["read_count", "umi_filtering", "adapter_trimmed", "joined"]
    messages = ["unprocessed fastqs", "umi filtered fastqs", "trimmed fastqs", "joined fastqs"]
    process_order = [unprocessed_fastqs, umi_filtered_fastqs, adapter_trimmed_fastqs, joined_fastqs]
    metadata = defaultdict(dict)

    for key, message, file_group in zip(out_order, messages, process_order):
        process_samples(file_group, message, metadata, key)

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
