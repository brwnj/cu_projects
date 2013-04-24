#!/usr/bin/env python
# encoding: utf-8

def main(args):
    fraction = "Fraction"
    for line in open(args.CSV, 'r'):
        line = line.strip().split(",")
        try:
            distance = line[0].split("(")[0]
            absorbance = line[1]
            try:
                float(line[2])
                fraction = line[2]
            except ValueError:
                pass
            print "\t".join([distance, absorbance, fraction])
        except IndexError:
            continue

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("CSV", help="gradient file")
    args = p.parse_args()
    main(args)