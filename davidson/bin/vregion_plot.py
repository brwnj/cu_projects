#!/usr/bin/env python
# encoding: utf-8
"""
"""
from itertools import groupby
from toolshed import nopen
from collections import Counter

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
    print "v_region\tcoverage"
    for v_region, count in sorted(v_regions.iteritems()):
        print "%s\t%s" % (v_region, count)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fasta", help="filtered fasta with detailed headers")
    main(p.parse_args())