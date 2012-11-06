#!/usr/bin/env python
# encoding: utf-8
"""
Trim peaks by finding the point of maximum intensity then expanding to the
size of width.

chr1	21877673	21877674	p.84196.3|USP48	1000	+
chr1	21877989	21877990	p.84196.2|USP48	1000	+
chr1	21878312	21878313	p.84196.1|USP48	1000	+
chr1	22188427	22188428	p.23436.1|CELA3B	1000	+
chr1	22211615	22211616	p.10136.1|CELA3A	1000	+
chr1	22211615	22211616	p.23436.3|CELA3A	1000	+
chr1	22289206	22289207	p.998.1|CDC42	1000	+
chr1	22289878	22289879	p.998.2|CDC42	1000	+
chr1	22290690	22290691	p.998.3|CDC42	1000	+
chr1	22291427	22291428	p.998.4|CDC42	1000	+
chr1	22291490	22291491	p.998.5|CDC42	1000	+
chr1	22291811	22291812	p.998.6|CDC42	1000	+
chr1	22292020	22292021	p.998.7|CDC42	1000	+

python ../../bin/cleavage_sites.py --annotation ~/ref/hg18/polyautr3intersection.bed --genomedata ~/projects/polya/data/20121010/gdarch --tracka MP51.umi --trackb MP52.umi

real site:
PAFAH1B1

"""

import sys
import numpy as np
import scipy.stats as ss
from genomedata import Genome
from pybedtools import BedTool
from toolshed import reader
from collections import defaultdict, OrderedDict
from itertools import izip

import cPickle as pickle
import os.path as op


def maxvals(bedtool, genome, atrack, btrack, cutoff, scalea, scaleb):
    # poly(a) site window defined as 20 bp
    start = bedtool.start - 20
    stop = bedtool.stop + 20
    
    a = genome[bedtool.chrom][start:stop, atrack]
    b = genome[bedtool.chrom][start:stop, btrack]
    
    maxa = np.nanmax(a) * scalea
    maxb = np.nanmax(b) * scaleb
    
    # normalizing peaks with few or no reads
    if np.isnan(maxa):# or maxa < cutoff:
        maxa = 1.
    if np.isnan(maxb):# or maxb < cutoff:
        maxb = 1.
    
    return maxa, maxb


def scalingfactor(a, b):
    """returns scaling factor for each sample"""
    a = float(a)
    if a > b:
        return b / a, 1.0
    else:
        return 1.0, a / b


def main(args):
    annotation = BedTool(args.annotation)
    genome = Genome(args.genomedata)
    
    libsize = dict(zip(genome.tracknames_continuous, genome.sums))
    sfa, sfb = scalingfactor(libsize[args.tracka], libsize[args.trackb])
    
    autrs = defaultdict(dict)
    butrs = defaultdict(dict)
    
    # each line of the annotation is a poly(a) site
    for i, line in enumerate(annotation):
        # if i > 1000: break
        # parallelize by chrom if supplied
        if args.chrom and line.chrom != args.chrom: continue
        pola, polb = maxvals(line, genome, args.tracka, args.trackb, args.cutoff, sfa, sfb)
        polyaname, gene = line.name.split("|")
        autrs[gene].update({polyaname:pola})
        butrs[gene].update({polyaname:polb})
        
    # autrs['test'].update({"test1":15.})
    # autrs['test'].update({"test2":20.})
    # butrs['test'].update({"test1":25.})
    # butrs['test'].update({"test2":15.})
    
    for gene in autrs.keys():
        
        # dict of poly(a)
        pola = autrs.get(gene)
        polb = butrs.get(gene)
        
        # order poly(a) sites by decreasing peak heights
        pola = OrderedDict(sorted(pola.items(), key=lambda t: t[1]))
        polb = OrderedDict(sorted(polb.items(), key=lambda t: t[1]))
        
        # loop over the keys
        for aname, bname in izip(pola, polb):
            if aname == bname: continue
            if pola.get(aname) < args.cutoff and polb.get(bname) < args.cutoff: continue
            polastr = ""
            polbstr = ""
            for ka, va in pola.iteritems():
                polastr = polastr + ka + "\t" + str(va) + "\t"
            for kb, vb in polb.iteritems():
                polbstr = polbstr + kb + "\t" + str(vb) + "\t"
            
            print gene
            print polastr.rstrip("\t")
            print polbstr.rstrip("\t")
            break


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotation', required=True, help="3' UTR poly(A) sites as bed")
    p.add_argument('--genomedata', required=True, help="genome data archive")
    p.add_argument('--tracka', required=True, help="track name a")
    p.add_argument('--trackb', required=True, help="track name b")
    p.add_argument('--chrom', action='store_true', default=False, help='specify to split jobs based on chrom name')
    p.add_argument('--cutoff', default=7, type=int, help="minimum allowable peak height")
    args = p.parse_args()
    main(args)