#!/usr/bin/env python
# encoding: utf-8
"""
Create new fastq of the unmapped reads.
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

def dict_of_mapped(fp):
    s = set()
    for name in nopen(fp):
        s.add(name[:-1])
    return s

def main(args):
    mapped = dict_of_mapped(args.txt)
    with nopen(args.fastq) as fh:
        for name, seq, qual in readfx(fh):
            if name in mapped: continue
            print "@%s\n%s\n+\n%s" % (name, seq, qual)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("fastq", help="trimmed reads")
    p.add_argument("txt", help="mapped reads by name only")
    args = p.parse_args()
    main(args)
