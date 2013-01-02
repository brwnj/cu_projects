#!/usr/bin/env python
# encoding: utf-8
"""
Join R1 and R2 into SSAKE compatible fasta.
"""
import sys
from toolshed import nopen

class FastqReader(object):
    """Yields name, seq, qual."""
    def __init__(self, fastq):
        self.fastq = nopen(fastq)

    def __iter__(self):
        fq = self.fastq
        while True:
            id1 = fq.next().strip()
            seq = fq.next().strip()
            id2 = fq.next().strip()
            qual = fq.next().strip()
            if qual == "":
                if id1 != "":
                    sys.stderr.write(">> Incomplete fastq... skipping.\n")
                break
            yield id1[1:], seq, qual

def fastqtodict(fastq, verbose):
    if verbose:
        sys.stderr.write(">> Creating R2 reference...\n")
    fdict = {}
    for name, seq, qual in FastqReader(fastq):
        fdict[name] = seq
    return fdict

def main(args):
    r2 = fastqtodict(args.R2, args.verbose)
    if args.verbose:
        sys.stderr.write(">> Joining reads...\n")
    for name, seq, qual in FastqReader(args.R1):
        try:
            r2seq = r2.get(name)
            print ">%s:%d\n%s:%s" % (name, args.insert, seq, r2seq)
        except KeyError:
            sys.stderr.write(">> No match found for: %s\n" % read)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('R1')
    p.add_argument('R2')
    p.add_argument('--insert', type=int, default = 200,
            help="target insert size")
    p.add_argument('--verbose', '-v', action="store_true", 
            help="maximum verbosity")
    main(p.parse_args())