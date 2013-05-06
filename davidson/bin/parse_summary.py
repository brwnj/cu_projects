#!/usr/bin/env python
# encoding: utf-8
"""
"""
import sys
from toolshed import reader
from collections import Counter, defaultdict

def get_unique_protein_seqs(imgt):
    """docstring for get_unique_protein_seqs"""
    c = Counter()
    for l in reader(imgt, header=True):
        if "productive" not in l['Functionality']: continue
        if len(l['AA JUNCTION']) < 2: continue
        c.update([l['AA JUNCTION']])
    return c

def get_vdj_regions(counter, imgt):
    """docstring for get_vdj_regions"""
    p = defaultdict(list)
    s = {}
    for l in reader(imgt, header=True):
        if "productive" not in l['Functionality']: continue
        try:
            v = l["V-GENE and allele"].split()[1]
        except IndexError:
            v = "na"
        try:
            j = l["J-GENE and allele"].split()[1]
        except IndexError:
            j = "na"
        try:
            d = l["D-GENE and allele"].split()[1]
        except IndexError:
            d = "na"
        composition = "%s,%s,%s" % (v, d, j)
        protein_seq = l["AA JUNCTION"]
        p[protein_seq].append(composition)
        try:
            if len(l['Sequence']) > len(s[protein_seq]):
                s[protein_seq] = l['Sequence']
        except KeyError:
            s[protein_seq] = l['Sequence']
    return p, s

def main(args):
    seqs = get_unique_protein_seqs(args.imgt)
    vdj, longest_seqs = get_vdj_regions(seqs, args.imgt)
    for (seq, count) in seqs.most_common():
        fields = [seq, count]
        vdj_info = ""
        c = Counter()
        for combo in vdj[seq]:
            c.update([combo])
        for (combo, count) in c.most_common():
            vdj_info += "%s (%d);" % (combo, count)
        fields.append(vdj_info)
        fields.append(longest_seqs[seq].upper())
        print "\t".join(map(str, fields))

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("imgt", help="imgt report summary")
    args = p.parse_args()
    main(args)
