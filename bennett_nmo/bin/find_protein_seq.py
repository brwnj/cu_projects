#!/usr/bin/env python
# encoding: utf-8
"""
YYC, YFC, and YSC
"""
import sys
from Bio import trie
# from Bio.Seq import Seq
from toolshed import reader
from itertools import ifilterfalse

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
            if not l['CDR3-IMGT'].strip(): continue
            seqs[l['Sequence ID']] = l['CDR3-IMGT']
        except KeyError:
            continue
    return seqs

def chunker(it, n):
    # chunker('AAAABBBC', 4) --> AAAA AAAB AABB ABBB BBBC
    return [it[i:i+n] for i in xrange(0, len(it)+1-n, 1)]

def build_trie(target, query):
    t = trie.trie()
    lengths = sorted(set([len(v) for v in target.values()]))
    for name, seq in query.iteritems():
        seq_len = len(seq)
        for l in lengths:
            if l > seq_len: continue
            for subseq in chunker(seq, l):
                if t.has_key(subseq):
                    t[subseq] += ";%s" % name
                t[subseq] = name
    return t

def unique_everseen(iterable, key=None):
    seen = set()
    seen_add = seen.add
    if key is None:
        for element in ifilterfalse(seen.__contains__, iterable):
            seen_add(element)
            yield element
    else:
        for element in iterable:
            k = key(element)
            if k not in seen:
                seen_add(k)
                yield element

def process_seqs(tree, query):
    for name, seq in unique_everseen(query.iteritems(), key=lambda t: t[1]):
        if len(seq) < 5: continue
        print seq
        print tree.get_approximate(seq, 2)

def main(cdr3, imgtaa):
    proteins = get_proteins(cdr3)
    seqs = get_seqs(imgtaa)
    t = build_trie(proteins, seqs)
    process_seqs(t, proteins)
    
    # check if found cdr3 are subseqs of "known"
    t = build_trie(seqs, proteins)
    process_seqs(t, seqs)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cdr3", help="id and protein sequence file")
    p.add_argument("imgtaa", help="processed, ungapped aa via high-v")
    args = vars(p.parse_args())
    main(**args)
