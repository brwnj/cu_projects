#!/usr/bin/env python
# encoding: utf-8
"""
Combine IMGT AA sequences results into single files for each sample.
"""
import sys
import argparse
import os.path as op
from toolshed import reader, header

def main(result_files):
    common_header = header(result_files[0])[:-1]
    samples = set()
    for file in result_files:
        samples.add(op.basename(file).rsplit("_", 2)[0])
    # leaves one entry per unique sample
    for sample in samples:
        # open a new file to facilitate the join
        out = open("{sample}.combined.txt".format(sample=sample), "wb")
        i = 1
        print >>out, "\t".join(common_header)
        for result_file in result_files:
            if not op.basename(result_file).startswith(sample): continue
            for toks in reader(result_file, header=common_header):
                if not toks['Functionality'] == 'productive': continue
                toks['Sequence number'] = str(i)
                print >>out, "\t".join(toks[h] for h in common_header)
                i += 1
        out.close()

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("RESULT_FILES", nargs="+", help="IMGT AA result files")
    args = p.parse_args()
    main(args.RESULT_FILES)
