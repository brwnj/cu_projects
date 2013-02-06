#!/usr/bin/env python
# encoding: utf-8
import sys
import pandas as pd
from toolshed import nopen
import matplotlib.pyplot as plt

def main(args):
    for i, line in enumerate(nopen(args.CSV)):
        line = line.split(",")
        try:
            distance = float(line[0])
            absorbance = float(line[1])
            if line[2] and args.fraction:
                if distance == 0:
                    continue
                t = float(line[2])
                print "%f\t0" % (distance)
            else:
                print "%f\t%f" % (distance, absorbance)
        except ValueError:
            continue

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("CSV", help="gradient file")
    p.add_argument("-f", "--fraction", action="store_true")
    main(p.parse_args())