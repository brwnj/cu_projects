#!/usr/bin/env python
# encoding: utf-8
"""
Compare any column of two files.
TODO: fix bonly

ars_bl6 [ 857 [ 422 ] 321 ]

"""
import argparse
import sys
from toolshed import reader, header

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def get_args():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('filea')
    p.add_argument('fileb')
    p.add_argument('-m', '--mode', required=True,
                    choices="intersection union aonly bonly".split(),
                    default='intersection', help='type of output')
    p.add_argument('-col', '--column', required=True, 
                    help='column name if header, number if no header')
    p.add_argument('-sep', '--separator', default="\t", 
                    help='field delimiter')
    args = p.parse_args()
    return args


def get_set(filein, column, sep="\t"):
    fset = []
    if column.isdigit():
        column = int(column) - 1
        for line in reader(filein, header=False, sep=sep):
            fset.append(line[column])
    else:
        for line in reader(filein, header=True, sep=sep):
            fset.append(line[column])
    return set(fset)


def get_dict(filein, column, sep="\t"):
    fdict = {}
    if column.isdigit():
        column = int(column) - 1
        for line in reader(filein, header=False, sep=sep):
            fdict[line[column]] = line
    else:
        for line in reader(filein, header=True, sep=sep):
            fdict[line[column]] = line
    return fdict
    

def main():
    args = get_args()
    
    seta = get_set(args.filea, args.column, args.separator)
    setb = get_set(args.fileb, args.column, args.separator)
    
    comparisons = {'aonly':seta - setb,
                   'bonly':setb - seta,
                   'intersection': seta & setb,
                   'union':seta | setb}
    
    lookup = get_dict(args.filea, args.column, args.separator)
    
    if args.column.isdigit():
        for setitem in comparisons[args.mode]:
            print "\t".join([i for i in lookup[setitem]])
    else:
        head = header(args.filea)
        for setitem in comparisons[args.mode]:
            print "\t".join([lookup[setitem][i] for i in head])


if __name__ == "__main__":
    main()