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

"""

import sys
import numpy as np
import scipy.stats as ss
from genomedata import Genome
from pybedtools import BedTool
from toolshed import reader
from collections import defaultdict


def trackvals(bedtool, genome, controltrack, treattrack, sfactora, sfactorb):
    start = bedtool.start - 10
    stop = bedtool.stop + 10
    
    control = genome[bedtool.chrom][start:stop, controltrack]
    treat = genome[bedtool.chrom][start:stop, treattrack]
    
    scaledcontrol = np.nanmax(control) * sfactora
    scaledtreat = np.nanmax(treat) * sfactorb
    
    # if scaledcontrol < cutoff and scaledtreat < cutoff:
    #     return 1.0
    # 
    # p = ss.binom_test([scaledcontrol, scaledtreat])
    # 
    # return p if p > 0.0 else 1.0
    
    return scaledcontrol, scaledtreat

#normalizing a matrix
df_norm = (df - df.mean()) / (df.max() - df.min())

def main(args):
    annotation = BedTool(args.annotation)
    genome = Genome(args.genomedata)
    
    # get the values into a dataframe - the future will have replicates
    # normalize the data in the dataframe
    # binomial test and resultant pvalue for each known polya site
    # with the pvalues, run qvality
    # loop through 3'UTRs, saving only ones having multiple with significant qvalues 
    
    # pvals = defaultdict(list)
    avals = defaultdict(list)
    bvals = defaultdict(list)
    
    for line in annotation:
        # parallelize by chrom if supplied
        if args.chrom and line.chrom != args.chrom: continue
        # get p-value
        # p = binomialtest(b, genome, args.tracka, args.trackb, sfa, sfb, args.cutoff)
        a, b = trackvals(line, genome, args.tracka, args.trackb, sfa, sfb)
        gene = line.name.split("|")[1]
        # pvals[gene].append(p)
        # avals[gene].append(a)
        bvals[gene].append(b)
    
    for gene in avals.keys():



    # for gene, results in pvals.iteritems():
        # pvals < 0.05
        # sig = [v for v in results if 0.0 < v < 0.05]
        # print gene name, all pvals, count of sig pvals
        # print gene + "\t" + ",".join(map(str, results)) + "\t" + str(len(sig))

        # obs = np.array(avals.get(gene))
        # exp = np.array(bvals.get(gene)) * np.sum(obs)
        # p = ss.chisquare(obs, exp)[1]
        # print gene + "\t" + str(p)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--annotation', required=True, help="3' UTR poly(A) sites as bed")
    p.add_argument('--genomedata', required=True, help="genome data archive")
    p.add_argument('--tracka', required=True, help="track name a")
    p.add_argument('--trackb', required=True, help="track name b")
    p.add_argument('--chrom', action='store_true', default=False, help='specify to split jobs based on chrom name')
    p.add_argument('--cutoff', default=5, type=int, help="minimum allowable peak height")
    args = p.parse_args()
    main(args)