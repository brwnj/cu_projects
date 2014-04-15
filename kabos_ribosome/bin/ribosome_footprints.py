#!/usr/bin/env python
# coding=utf-8
"""
"""
from toolshed import reader
from collections import Counter
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter


class Bed(object):
    __slots__ = ['chr','start','stop','name','score','strand']
    def __init__(self, args):
        for k, v in zip(self.__slots__, args):
            setattr(self, k, v)
        self.start = int(self.start)
        self.stop = int(self.stop)
        self.score = int(self.score)

    def rev(self):
        return True if self.strand == "-" else False


class BedGraph(object):
    def __init__(self, args):
        self.chr = args[0]
        self.start = args[1]
        self.stop = args[2]
        self.count = int(float(args[3]))


def bedgraph_to_dict(bedgraph):
    # d[chr":"start":"stop] = count
    d = {}
    for bg in reader(bedgraph, header=BedGraph):
        d["%s:%s:%s" % (bg.chr, bg.start, bg.stop)] = bg.count
    return d


def main(ribosome_bg, mrna_bg, startcodon_bed):
    codon_span = 12

    counts = {}
    counts['rrna'] = Counter()
    counts['mrna'] = Counter()

    rrna = bedgraph_to_dict(ribosome_bg)
    mrna = bedgraph_to_dict(mrna_bg)

    for b in reader(startcodon_bed, header=Bed):
        for i in xrange(0, codon_span):
            mult = i * 3

            if b.rev:
                zero = "%s:%d:%d" % (b.chr, b.stop - 1 - mult, b.stop - 0 - mult)
                one = "%s:%d:%d" % (b.chr, b.stop - 2 - mult, b.stop - 1 - mult)
                two = "%s:%d:%d" % (b.chr, b.stop - 3 - mult, b.stop - 2 - mult)

            else:
                zero = "%s:%d:%d" % (b.chr, b.start + 0 + mult, b.start + 1 + mult)
                one = "%s:%d:%d" % (b.chr, b.start + 1 + mult, b.start + 2 + mult)
                two = "%s:%d:%d" % (b.chr, b.start + 2 + mult, b.start + 3 + mult)

            try:
                counts['rrna']['zero'] += rrna[zero]
            except KeyError:
                pass
            try:
                counts['mrna']['zero'] += mrna[zero]
            except KeyError:
                pass
            try:
                counts['rrna']['one'] += rrna[one]
            except KeyError:
                pass
            try:
                counts['mrna']['one'] += mrna[one]
            except KeyError:
                pass
            try:
                counts['rrna']['two'] += rrna[two]
            except KeyError:
                pass
            try:
                counts['mrna']['two'] += mrna[two]
            except KeyError:
                pass

    print "\t".join(["pos", "rrna", "mrna"])
    print "0\t%d\t%d" % (counts['rrna']['zero'], counts['mrna']['zero'])
    print "1\t%d\t%d" % (counts['rrna']['one'], counts['mrna']['one'])
    print "2\t%d\t%d" % (counts['rrna']['two'], counts['mrna']['two'])


if __name__ == '__main__':
    p = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('ribosome_bg', help="ribosome fraction bedgraph")
    p.add_argument('mrna_bg', help="mRNA bedgraph")
    p.add_argument('startcodon_bed', help="unique start codon sites as bed")
    args = vars(p.parse_args())
    main(**args)
