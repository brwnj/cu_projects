#!/usr/bin/env python
# encoding: utf-8
"""
convert scarf to fastq.
"""
from toolshed import nopen

def main(args):
    for scarfline in nopen(args.scarf):
        name, seq, qual = scarfline.strip().rsplit(":", 2)
        assert(len(seq) == len(qual))
        print "@%s\n%s\n+\n%s" % (name, seq, qual)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('scarf')
    args = p.parse_args()
    main(args)