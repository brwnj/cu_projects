#!/usr/bin/env python
# encoding: utf-8
"""
for f in *; do
    echo $f;
    echo "bedtools intersect -f .25 -c -a ~/projects/hits-clip/data/common/mirbase/mirbase-18/hsa.hg18.bed.gz -b ${f}/${f}.bed.gz | gzip -c > ${f}/${f}.mirna_abundance.bed.gz" \
        | bsub-ez ${f}.btinters -q idle
done
"""
import sys
import os.path as op
from toolshed import reader
from collections import defaultdict

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def main(args):
    gm = defaultdict(dict)
    for f in args.files:
        for l in reader(f, header="chrom start stop name score strand abundance".split()):
            try:
                if int(l['abundance']) == 0: continue
                # gm[<parsed file name>][<miRNA name>] = abundance value
                gm[op.basename(f).split(".", 1)[0]][l['name']] = l['abundance']
            except KeyError:
                # header failed to set l['abundance']
                pass

    # the sample names
    caselist = sorted(gm.keys())
    # only save lines where at least one sample has a positive value
    completeset = []
    for i, case in enumerate(caselist):
        keys = sorted(gm[caselist[i]].keys())
        for k in keys:
            completeset.append(k)
    mirnas = set(completeset)
    
    # print the matrix
    print "#mirna_id\t" + "\t".join(k for k in caselist)
    for mirna in mirnas:
        fields = [mirna]
        for c in caselist:
            try:
                fields.append(gm[c][mirna])
            except KeyError:
                # miRNA not present in this case
                fields.append("0.0")
        print "\t".join(map(str, fields))


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('files', nargs='+', 
            help='space delimited list or glob')

    args = p.parse_args()
    main(args)