#!/usr/bin/env python
# encoding: utf-8
"""
"""
import os
import argparse

def gsamples(fname):
    fname = fname.split(".")[0].split("_")
    return int(fname[0][2:]), int(fname[2][2:])

def main(d):
    for f in os.listdir(d):
        a, b = gsamples(f)
        if a > b:
            os.remove(f)

if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('direct', help="place to fix")
    args = p.parse_args()
    main(args.direct)
