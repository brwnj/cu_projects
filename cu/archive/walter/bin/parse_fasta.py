#!/usr/bin/env python
"""
Parse UCSC downloaded fasta to include only transcript name and strand information
in the name line.

Original line:
>hg19_ensGene_ENST00000237247 range=chr1:67208779-67210057 5'pad=0 3'pad=0 strand=+ repeatMasking=none

Result:
>ENST00000237247_+
"""
import sys
import itertools
from toolshed import nopen

def read_fasta(fa):
    fh = nopen(fa)
    for header, group in itertools.groupby(fh, lambda line: line[0] == '>'):
        if header:
            line = group.next()
            info = line[1:].strip().split()
            transcript = info[0].rsplit("_", 1)[1]
            strand = info[4].split("=")[1]
            name = "%s_%s" % (transcript, strand)
        else:
            seq = ''.join(line.strip() for line in group)
            yield name, seq

def print_fasta(name, seq, wrap=70):
    print ">%s" % name
    # sequence lines standardized at length of wrap
    print "\n".join([seq[i:i + wrap] for i in range(0, len(seq), wrap)])

def main(args):
    for name, seq in read_fasta(args.fasta):
        print_fasta(name, seq)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('fasta', help='UCSC downloaded fasta')
    main(p.parse_args())