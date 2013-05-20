#!/usr/bin/env python
# encoding: utf-8
"""
Add piece or all of index read back onto R1 (or R2).
"""
from itertools import izip
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
        name, seqs, last = last[1:], [], None
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

def main(args):
    with nopen(args.index) as idx, nopen(args.fastq) as fq:
        for (xname, xseq, xqual), (name, seq, qual) in izip(readfx(idx), readfx(fq)):
            xname = xname.partition(" ")[0]
            xref = name.partition(" ")[0]

            assert xname == xref
            
            if args.end:
                xseq = xseq[args.begin:args.end]
                xqual = xqual[args.begin:args.end]
            else:
                xseq = xseq[args.begin:]
                xqual = xqual[args.begin:]
            
            print "@{name}\n{xseq}{seq}\n+\n{xqual}{qual}".format(name=name,
                            xseq=xseq, seq=seq, xqual=xqual, qual=qual)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("index", metavar="INDEX", help="index read")
    p.add_argument("fastq", metavar="FASTQ", help="R1 or R2")
    p.add_argument("-b", dest="begin", default=0, type=int,
            help="0-based start of index read to save [%(default)s]")
    p.add_argument("-e", dest="end", type=int,
            help="0-based end of index read to save")
    args = p.parse_args()
    main(args)
