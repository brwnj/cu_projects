#!/usr/bin/env python
# encoding: utf-8
r"""convert deseq output (csv) to bed format

python result2bed result_t1.csv | bedtools sort -i - > result_t1.beds
"""
import sys
import tempfile
from toolshed import reader
from collections import defaultdict

def main(args):
    for r in reader(args.result, header=True, sep=","):
        if r['pval'] == "NA" or float(r['pval']) > args.pvalue: continue
        chrom, coords = r['id'].split(':')
        start, stop = coords.split('-')
        fields = [chrom, start, stop, r['id'], ".", r['pval']]
        # not sorted
        print "\t".join(map(str, fields))

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('result', help='DESeq result file')
    p.add_argument('-p','--pvalue', default=0.05, type=float, help='desire padj cutoff value [ 0.05 ]')
    args = p.parse_args()
    main(args)