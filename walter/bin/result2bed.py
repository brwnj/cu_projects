#!/usr/bin/env python
# encoding: utf-8
"""
i messed up the input and didn't put coords in the miRNA name field.

this is the fix.
"""
import sys
import tempfile
from toolshed import reader
from collections import defaultdict

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def buildpeakcoords(bed):
    """build a simple lookup dict since i lost coord info in deseq"""
    pc = defaultdict(dict)
    for b in reader(bed, header="chrom start stop name".split()):
        pc[b['name']]['chrom'] = b['chrom']
        pc[b['name']]['start'] = b['start']
        pc[b['name']]['stop'] = b['stop']
    return pc


def main(args):
    # build dictionary of consensus peak information
    coords = buildpeakcoords(args.consensus)
    # convert results to bed saving only padj in score field
    f = open(tempfile.mktemp(suffix=".deseq"), "w")
    for r in reader(args.result, header=True, sep=","):
        #if r['padj'] == "NA" or float(r['padj']) > args.pvalue: continue
        peak = coords.get(r['id'])
        fields = [peak.get('chrom'), peak.get('start'), peak.get('stop'), r['id'], r['padj']]
        f.write("\t".join(map(str, fields)) + "\n")
    f.close()
    for t in reader(f.name, header=False):
        print "\t".join(map(str, t))


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('result', help='DESeq result file')
    p.add_argument('consensus', help='consensus peak bed')
    p.add_argument('gff', help='mirbase annotation as gff')
    p.add_argument('-p','--pvalue', default=0.01, type=float, help='desire padj cutoff value [ 0.01 ]')
    args = p.parse_args()
    main(args)