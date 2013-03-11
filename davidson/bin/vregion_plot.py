#!/usr/bin/env python
# encoding: utf-8
"""
"""
from itertools import groupby
from toolshed import nopen
from collections import Counter

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

def read_fasta(fa):
    with nopen(fa) as fh:
        for header, group in groupby(fh, lambda line: line[0] == '>'):
            if header:
                line = group.next()
                name = line[1:].strip()
            else:
                seq = ''.join(line.strip() for line in group)
                yield name, seq

def main(args):
    name_key = "contig length reads coverage seed vregion jregion".split()
    v_regions = Counter()
    with nopen(args.fasta) as fasta:
        for name, seq in read_fasta(fasta):
            # remove some text from iSSAKE output
            name = name.replace("size","").replace("cov","").replace("read","").replace("seed:","")
            name_d = dict(zip(name_key, name.split("|")))
            v_regions.update([name_d['vregion']])
    # stick a header on there
    # print "v_region\tcoverage"
    coverage = []
    regions = []
    for v_region, count in sorted(v_regions.iteritems()):
        # print "%s\t%s" % (v_region, count)
        regions.append(v_region)
        coverage.append(count)

    fig = plt.figure(figsize=(15,10), dpi=300)
    ax = fig.add_subplot(111)
    N = np.arange(len(coverage))
    bar_width = 2
    ax.bar(N, coverage, facecolor='red')
    ax.set_ylabel('Count')
    plt.xticks(rotation=70)
    ax.set_xticks(N)
    ax.set_xticklabels(regions)
    
    plt.title(args.title)
    plt.savefig(args.out, bbox_inches='tight', dpi=300)
    

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta", help="filtered fasta with detailed headers")
    p.add_argument("out", help="save plot as...")
    p.add_argument("--title", help="plot title")
    main(p.parse_args())