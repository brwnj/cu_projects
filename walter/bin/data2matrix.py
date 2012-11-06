#!/usr/bin/env python
# encoding: utf-8
"""
"""
import sys
from toolshed import reader
from collections import defaultdict

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def main(args):
    gm = defaultdict(dict)
    for f in args.files:
        for l in reader(f, header="chrom start stop name counts nonzero blength nonzerofracofb".split()):
            try:
                # gm[<parsed file name>][<peak name>] = count value
                fullname = "%s:%s:%s:%s" % (l['name'], l['chrom'], l['start'], l['stop'])
                gm[f.split(".", 1)[0]][fullname] = l['counts']
            except KeyError:
                # header failed to set l['val']
                pass

    # print the matrix
    caselist = sorted(gm.keys())
    # this step is unnecessary as they all have counts for the same peaks
    completeset = []
    for i, case in enumerate(caselist):
        keys = sorted(gm[caselist[i]].keys())
        for k in keys:
            completeset.append(k)
    
    peaks = set(completeset)
    print "#peak_name\t" + "\t".join(k for k in caselist)
    for peak in peaks:
        fields = [peak]
        for c in caselist:
            try:
                fields.append(gm[c][peak])
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