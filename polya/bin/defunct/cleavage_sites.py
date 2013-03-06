#!/usr/bin/env python
# encoding: utf-8
"""
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
import pandas as pd
import scipy.stats as ss
from genomedata import Genome
from pybedtools import BedTool
from collections import defaultdict

def getmax(bedtool, genome, track_name):
    # poly(a) site window defined as 20 bp
    start = bedtool.start - 20
    stop = bedtool.stop + 20
    
    maxval = np.nanmax(genome[bedtool.chrom][start:stop, track_name])
    if np.isnan(maxval):
        maxval = 0.
    return maxval

def main(args):
    annotation = BedTool(args.annotation)
    genome = Genome(args.genomedata)
    tracks = args.tracks.split(",")
    
    # container of peak values per polya per track
    maxvals = {}
    # list of poly(a)s per gene
    genes = defaultdict(list)
    for track in tracks:
        maxvals[track] = {}
        # each line of the annotation is a poly(a) site
        for bt in BedTool(args.annotation):
            # parallelize by chrom if supplied
            if args.chrom and bt.chrom != args.chrom: continue
            # format outlined in help
            polya_name, gene_name = bt.name.split("|")
            maxvals[track][polya_name] = getmax(bt, genome, track)
            genes[gene_name].append(polya_name)
    # dataframe of polya and peak heights
    df = pd.DataFrame(maxvals)
    
    # print df.to_string()
    
    # sys.exit(1)
    # normalize the data
    # df_norm = (df - df.mean()) / (df.max() - df.min())
    # print df_norm.to_string()
    pvals = defaultdict(dict)
    for polya, values in df.iterrows():
        a = np.array([values['MP51.umi']])
        b = np.array([values['MP52.umi']])
        p = ss.binom_test([a,b])
        pvals['p'].update({polya:p})
    pvals_df = pd.DataFrame(pvals)
    # print pvals_df.to_string()
    df = df.join(pvals_df)
    print df.to_string()
        
# 
# "p.93974.1"     44      36
# "p.9403.1"      27      20
# "p.9403.2"      0       1
# "p.9406.2"      2       0
# "p.9406.3"      0       0
# "p.9406.4"      108     142
# "p.9406.5"      108     142
# "p.9410.1"      14      24
# "p.9410.2"      0       0
# "p.94161.3"     2604    2738

    # binomial test and resultant pvalue for each known polya site
    # with the pvalues, run qvality
    # loop through 3'UTRs, saving only ones having multiple with significant qvalues 
    
    #binomial test = p = ss.binom_test([scaledcontrol, scaledtreat])
    
    # autrs['test'].update({"test1":15.})
    # autrs['test'].update({"test2":20.})
    # butrs['test'].update({"test1":25.})
    # butrs['test'].update({"test2":15.})
    
    # for gene in autrs.keys():
    #     
    #     # dict of poly(a)
    #     pola = autrs.get(gene)
    #     polb = butrs.get(gene)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('annotation',
            help="3' UTR poly(A) sites (bed); name column format is 'polya_name|gene_name'")
    p.add_argument('genomedata', help="genome data archive")
    p.add_argument('tracks', help="comma separated list of track names")
    p.add_argument('--chrom', action='store_true',
            help='specify to split jobs based on chrom name')
    args = p.parse_args()
    main(args)