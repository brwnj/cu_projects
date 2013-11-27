#!/usr/bin/env python
# encoding: utf-8
"""
fastq to IMGT compatible fastas.
"""

import os.path as op
from toolshed import nopen
from itertools import islice

def readfq(fq):
    class Fastq(object):
        def __init__(self, args):
            self.name = args[0][1:]
            self.seq = args[1]
            self.qual = args[3]
            assert len(self.seq) == len(self.qual)

        def __repr__(self):
            return "Fastq({name})".format(name=self.name)

        def __str__(self):
            return "@{name}\n{seq}\n+\n{qual}".format(name=self.name,
                    seq=self.seq, qual=self.qual)

    with nopen(fq) as fh:
        fqclean = (x.strip("\r\n") for x in fh if x.strip())
        while True:
            rd = [x for x in islice(fqclean, 4)]
            if not rd: raise StopIteration
            assert all(rd) and len(rd) == 4
            yield Fastq(rd)

def gfname(name, seg):
    return "{name}.{seg}.fasta".format(**locals())

def main(fastqs, n):
    out = ""
    for fastq in fastqs:
        seg = 1
        name = fastq.split(".f")[0]
        out = open(gfname(name, seg), 'wb')
        for i, r in enumerate(readfq(fastq), start=1):
            if i % n == 0:
                out.close()
                seg += 1
                out = open(gfname(name, seg), 'wb')
            out.write(">{name}\n{seq}\n".format(name=r.name, seq=r.seq))
        out.close()

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("fastqs", nargs="+", help="fastqs to convert")
    p.add_argument("-n", default=500000, type=int,
            help="number of lines per fasta")
    args = vars(p.parse_args())
    main(**args)
