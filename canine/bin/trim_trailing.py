#!/usr/bin/env python
# encoding: utf-8
"""
"""
from toolshed import nopen

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

def chunker(it, n):
    """
    >>> chunker('AAAABBBC', 4)
    ['AAAA', 'AAAB', 'AABB', 'ABBB', 'BBBC', 'BBC', 'BC', 'C']
    """
    return [it[i:i+n] for i in xrange(0, len(it), 1)]

def get_trim_loc(seq, base, fraction, wsize):
    for i, subseq in enumerate(chunker(seq, wsize)):
        if float(subseq.count(base))/len(subseq) >= fraction:
            return subseq.find(base) + i
    return len(seq)

def main(args):
    with nopen(args.fastq) as fh:
        for name, seq, qual in readfx(fh):
            pos = get_trim_loc(seq, args.b, args.f, args.w)
            if pos > args.l:
                print "@%s\n%s\n+\n%s" % (name, seq[:pos], qual[:pos])

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq")
    p.add_argument("-b", choices=['A','C','G','T'], required=True,
            help="base to trim")
    p.add_argument("-f", metavar="FRACTION", default=0.6, type=float,
            help="trimming threshold [%(default)s]")
    p.add_argument("-l", metavar="LENGTH", default=30, type=int,
            help="minimum sequence length [%(default)s]")
    p.add_argument("-w", metavar="WINDOW", default=15, type=int,
            help="window size [%(default)s]")
    args = p.parse_args()
    main(args)
