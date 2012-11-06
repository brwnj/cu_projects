#!/usr/bin/env python
# encoding: utf-8
"""
Compare sets of miRNAs of one group vs another. Either group may contain
multiple files.

Node input format:
hsa-miR-1976 = 0.0
hsa-miR-1973 = 54.0
hsa-miR-1972 = 0.0

Edge input format:
hsa-miR-502-5p (7) MAGT1 = 11.0
hsa-miR-5195-3p (7) MAGT1 = 11.0
hsa-miR-3680-3p (7) MAGT1 = 11.0

Node and Edge Outputs:
aname.bname.a_only.txt
aname.bname.b_only.txt
aname.bname.intersection.txt

sif input format:
hsa-miR-3183	8	RBMY1F
hsa-miR-1229	8	RBMY1F
hsa-miR-3190-5p	8	TTTY10

sif output:
HS5+HS27A vs hMSC
HS5+HS27A: 32227 nodes, sharing 3600 (11.17%) with hMSC
hMSC: 7466 nodes, sharing 3600 (48.22%) with HS5+HS27A
"""
import argparse
import sys
from toolshed import reader

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def get_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('-a', nargs='+', required=True, 
                    help='miRNA file(s) for group A')
    p.add_argument('-b', nargs='+', required=True,
                    help='miRNA file(s) for group B')
    p.add_argument('-aname', required=True,
                    help='name for group A')
    p.add_argument('-bname', required=True,
                    help='name for group B')
    p.add_argument('-mode', choices="node edge sif".split(), required=True,
                    help='type of file to process')
    args = p.parse_args()
    return args


def get_set_from_edge(filenamelist):
    """Takes file or list of files and returns a set of present miRNA seeds
    from a Cytoscape edge attributes file.
    """
    mirset = []
    for filename in filenamelist:
        for attr in reader(filename, \
                header="mirna length gene op intensity".split(), sep=" "):
            try:
                # any conditional statement would presumably need to be done
                # on normalized values which we are not currently using
                mirset.append("%s%s" % (attr['mirna'], attr['gene']))
            except KeyError:
                # file header
                pass
    return set(mirset)


def get_set_from_node(filenamelist):
    """Takes file or list of files and returns a set of present miRNAs from
    a Cytoscape node attributes file.
    """
    mirset = []
    for filename in filenamelist:
        for attr in reader(filename, header="mirna op value".split(), sep=" "):
            if attr['value'] == "0.0": continue
            mirset.append(attr['mirna'])
    return set(mirset)


def get_set_from_sif(sifs):
    """Takes file or list of files and returns a set of miRNAs from a 
    Cytoscape network file.
    """
    mirset = []
    for sif in sifs:
        for attr in reader(sif, header="mir length gene".split(), sep="\t"):
            mirset.append("%s|%s" % (attr['mir'], attr['gene']))
    return set(mirset)


def main():
    args = get_args()
    if args.mode is "node":
        seta = get_set_from_node(args.a)
        setb = get_set_from_node(args.b)
    elif args.mode is "edge":
        seta = get_set_from_edge(args.a)
        setb = get_set_from_edge(args.b)
    else:
        seta = get_set_from_sif(args.a)
        setb = get_set_from_sif(args.b)
    
    comparisons = {'%s_only' % args.aname:seta - setb,
                   '%s_only' % args.bname:setb - seta,
                   'intersection': seta & setb}
    
    if args.mode is "node" or args.mode is "edge":
        for test, result in comparisons.iteritems():
            fout = open("%s.%s.%s.txt" % (args.aname, args.bname, test), 'w')
            fout.write("\n".join(r for r in result))
    else:
        a_tot = len(seta)
        b_tot = len(setb)
        intersection = len(comparisons['intersection'])
        a_perc = intersection / float(a_tot) * 100
        b_perc = intersection / float(b_tot) * 100
        print "%s vs %s" % (args.aname, args.bname)
        print "%s: %d connections, sharing %d (%.2f%%) with %s" % \
            (args.aname, a_tot, intersection, a_perc, args.bname)
        print "%s: %d connections, sharing %d (%.2f%%) with %s" % \
            (args.bname, b_tot, intersection, b_perc, args.aname)


if __name__ == "__main__":
    main()
