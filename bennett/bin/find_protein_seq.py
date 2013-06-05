#!/usr/bin/env python
# encoding: utf-8
"""
YYC, YFC, and YSC
"""
import sys
# from Bio.Seq import Seq
from toolshed import reader
from Bio import trie

def get_proteins(f):
    d = {}
    for t in reader(f, header=["id", "pseq"]):
        # duplicates will be discarded
        d[t['id']] = t['pseq']
    return d

def get_seqs(fn):
    seqs = {}
    for l in reader(fn):
        try:
            if not l['AA JUNCTION'].strip(): continue
            seqs[l['Sequence ID']] = l['AA JUNCTION']
        except KeyError:
            continue
    return seqs

def chunker(it, n):
    # chunker('AAAABBBC', 4) --> AAAA AAAB AABB ABBB BBBC
    return [it[i:i+n] for i in xrange(0, len(it)+1-n, 1)]

def build_trie(p, seqs):
    t = trie.trie()
    lengths = sorted(set([len(v) for v in p.values()]))
    for name, seq in seqs.iteritems():
        seq_len = len(seq)
        for l in lengths:
            if l > seq_len: continue
            for subseq in chunker(seq, l):
                if t.has_key(subseq):
                    t[subseq] += ";%s" % name
                t[subseq] = name
    return t

def process_seqs(t, pseqs):
    for name, pseq in pseqs.iteritems():
        print pseq
        print t.get_approximate(pseq, 2)

def main(args):
    proteins = get_proteins(args.cdr3)
    seqs = get_seqs(args.seqs)
    t = build_trie(proteins, seqs)
    process_seqs(t, proteins)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cdr3", help="id and protein sequence file")
    p.add_argument("seqs", help="joined sequences")
    args = p.parse_args()
    main(args)
