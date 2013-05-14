#!/usr/bin/env python
# encoding: utf-8
"""
YYC, YFC, and YSC
"""
import sys
from Bio import trie, triefind
from Bio.Seq import Seq
from toolshed import reader, nopen

def readfx(fh):
    # https://github.com/lh3/readfq/blob/master/readfq.py
    last = None
    while True:
        if not last:
            for l in fh:
                if l[0] in '>@':
                    last = l[:-1]
                    break
        if not last: break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fh:
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':
            yield name, ''.join(seqs), None
            if not last: break
        else:
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fh:
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):
                    last = None
                    yield name, seq, ''.join(seqs);
                    break
            if last:
                yield name, seq, None
                break

def get_proteins(f):
    """docstring for get_proteins"""
    d = {}
    for t in reader(f, header=["id", "pseq"]):
        # duplicates will be discarded
        d[t['id']] = t['pseq']
    return d

def translate(seq):
    """translates the sequence in each reading frame"""
    seq = Seq(seq)
    i = -1
    while i < 2:
        i += 1
        yield seq[i:].translate(to_stop=True).tostring()

def chunker(it, n):
    # chunker('AAAABBBC', 4) --> AAAA AAAB AABB ABBB BBBC
    return [it[i:i+n] for i in xrange(0, len(it)+1-n, 1)]

def build_trie(p, seqs):
    """docstring for process_seqs"""
    t = trie.trie()
    lengths = sorted(set([len(v) for v in p.values()]))
    with nopen(seqs) as fh:
        for name, seq, qual in readfx(fh):
            for pseq in translate(seq):
                if "YYC" in pseq or "YFC" in pseq or "YSC" in pseq:
                    pseq_len = len(pseq)
                    for l in lengths:
                        if l > pseq_len: continue
                        for subseq in chunker(pseq, l):
                            if t.has_key(subseq):
                                t[subseq] += ";%s" % pseq
                            t[subseq] = name
                    break
    return t

def process_seqs(t, pseqs):
    #iterate through the protein sequences and see if there are any hits
    for uid, pseq in pseqs.iteritems():
        print pseq
        print t.get_approximate(pseq, 2)

def main(args):
    proteins = get_proteins(args.cdr3)
    t = build_trie(proteins, args.seqs)
    process_seqs(t, proteins)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("cdr3", help="id and protein sequence file")
    p.add_argument("seqs", help="joined sequences")
    args = p.parse_args()
    main(args)
