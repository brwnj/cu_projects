#!/usr/bin/env python
# encoding: utf-8
"""
Compare any column of two files.
"""
import sys
import urllib
import urllib2
from toolshed import reader, header
from itertools import izip


__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def get_set(filein, column, sep="\t"):
    fset = []
    if column.isdigit():
        for line in reader(filein, header=False, sep=sep):
            fset.append(line[int(column) - 1])
    else:
        for line in reader(filein, header=True, sep=sep):
            fset.append(line[column])
    return set(fset)


def get_dict(filein, column, sep="\t"):
    fdict = {}
    if column.isdigit():
        for line in reader(filein, header=False, sep=sep):
            fdict[line[int(column) - 1]] = line
    else:
        for line in reader(filein, header=True, sep=sep):
            fdict[line[column]] = line
    return fdict


def venn_gchart(a, b, c=None, colors=None, labels=None, size='300x300'):
    """
    a, b, and c are sets of the column being compared.
    *colors* is a list of 3 hex colors
    *labels* is a list of 3 labels
    *outfn* is the output PNG you want to create.
    *size* is the size in pixels for the PNG
    """
    # The order of values is meaningful to the API, see
    # http://code.google.com/apis/chart/docs/gallery/venn_charts.html
    if c:
        vals = [len(a),
                len(b),
                len(c),
                len(a&b),
                len(a&c),
                len(b&c),
                len(a&b&c)]
        if labels:
            print "Total size %s: %d" % (labels[0], vals[0])
            print "Total size %s: %d" % (labels[1], vals[1])
            print "Total size %s: %d" % (labels[2], vals[2])
            print "Overlap %s %s: %d" % (labels[0], labels[1], vals[3])
            print "Overlap %s %s: %d" % (labels[0], labels[2], vals[4])
            print "Overlap %s %s: %d" % (labels[1], labels[2], vals[5])
            print "Overlap %s %s %s: %d" % (labels[0], labels[1], labels[2], vals[6])

    else:
        # insert 0 for size of 3rd circle.
        vals = [len(a), len(b), 0, len(a&b)]
        if labels:
            print "Total size %s: %d" % (labels[0], vals[0])
            print "Total size %s: %d" % (labels[1], vals[1])
            print "Overlap %s %s: %d" % (labels[0], labels[1], vals[3])
        labels = labels[:2]
    # API doesn't seem to like large numbers, so get fractions instead, then
    # join make a comma-separated list of values.
    mx = float(max(vals))
    vals = [i/mx for i in vals]
    valstr = ','.join(map(str,vals))

    data = {'cht':'v',
            'chs':size,
            'chd':'t:'+valstr}

    # Add the optional data, if specified
    if labels:
        data['chdl'] = '|'.join(labels)
    if colors:
        data['chco'] = ','.join(colors)
    return data


def gchart(data, outfn='out.png'):
    """
    Sends data to Google Chart API
    """
    data = urllib.urlencode(data)

    url = 'https://chart.googleapis.com/chart?'

    # Request and get the PNG
    req = urllib2.Request(url, data)
    print url + data
    response = urllib2.urlopen(req)
    f = open(outfn,'w')
    f.write(response.read())
    f.close()


def main(args):
    sets = []
    args.a = get_set(args.a, args.column, args.sep)
    args.b = get_set(args.b, args.column, args.sep)
    if args.c:
        args.c = get_set(args.c, args.column, args.sep)
    
    data = venn_gchart(a=args.a, b=args.b, c=args.c, colors=args.colors.split(","), labels=args.labels.split(","), size=args.size)
    gchart(data, outfn=args.out)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('-a', '--filea', dest='a', required=True, help="file for left-most circle")
    p.add_argument('-b', '--fileb', dest='b', required=True, help="file for right-most circle")
    p.add_argument('-c', '--filec', dest='c', default=None, help="file for bottom circle")
    p.add_argument('-l', '--labels', required=True, help='comma separated list of labels for a, b, and c')
    p.add_argument('-s', '--size', default='300x300', help="optional size of image")
    p.add_argument('-o', '--out', default='out.png', help="output file name")
    p.add_argument('--colors', help='optional comma-separated list of hex colors for circles a, b, and c',
                    default='00FF00,FF0000,0000FF')
    p.add_argument('-col', '--column', default="1",
                    help='column name if header, number if no header')
    p.add_argument('-sep', '--separator', dest='sep', default="\t", 
                    help='field delimiter')
    args = p.parse_args()
    main(args)