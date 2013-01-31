#!/usr/bin/env python
"""
from input bigwigs, create track strings without a color selection.
"""
import sys
import os.path as op

def main(args):
    for bw in args.bw:
        f_name = op.basename(bw)
        s_name = op.splitext(f_name)[0]
        track = "track type=bigWig name={0} description={0} \
bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/{1}/{2}/{3} \
maxHeightPixels=15:50:35 color= visibility=full".format(s_name, args.pi, args.date, f_name)
        print track

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, 
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('pi', help="PI folder name")
    p.add_argument('date', help="date folder name")
    p.add_argument('bw', metavar="BW", nargs="+", help="bigwigs")
    main(p.parse_args())