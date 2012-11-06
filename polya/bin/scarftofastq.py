#!/usr/bin/env python
# encoding: utf-8
"""
convert scarf to fastq.
"""
import sys
from toolshed import nopen

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def main(args):
    """do something amazing"""
    for scarfline in nopen(args.scarf):
        name, seq, qual = scarfline.strip().rsplit(":",2)
        print "@%s\n%s\n+\n%s" % (name, seq, qual)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(__doc__)
    p.add_argument('scarf')
    args = p.parse_args()
    main(args)