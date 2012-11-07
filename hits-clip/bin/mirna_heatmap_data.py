#!/usr/bin/env python
# encoding: utf-8
"""
generate the lists from the cases, then grab the required data from the samples.
"""
# import sys
import os.path as op
import pandas as pd
# import matplotlib.pyplot as plt
import numpy as np
from toolshed import reader

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def main(args):
    sampletocase = {'BT474herc':['PK23'],
                    'HELA':['helaa', 'helab'],
                    'MCF7':['PK11', 'PK31', 'PK51'],
                    'MCF7estr':['PK12', 'PK32'],
                    'MDA231':['PK24', 'PK42', 'PK54'],
                    'BT474':['PK21', 'PK41', 'PK52'],
                    'BT474estr':['PK22', 'PK53'],
                    'HS27A':['MP1', 'MP21', 'MP35'],
                    'HS5':['MP2', 'MP20', 'MP34'],
                    'hMSC':['MP36', 'MP43.ACTG', 'MP43.TCGA', 'MP44.ACTG', 'MP44.TCGA'],
                    'BMEC':['MP42.ACTG', 'MP45.ACTG', 'MP45.TCGA'],
                    'HUVEC':['MP24', 'MP38']}

    # create a dataframe with all miRNA counts from all samples
    allcounts = {}
    for f in args.counts:
        fname = op.basename(f).split(args.fullext)[0]
        casecounts = {}
        for line in reader(f, header="chrom start stop name score strand count".split()):
            casecounts[line['name']] = int(line['count'])
        allcounts[fname] = casecounts
    countsdf = pd.DataFrame(allcounts)

    # create a set of unique miRNAs from all the miRNA lists
    uniquemirnas = []
    for f in args.mirnalist:
        for line in reader(f, header=['name']):
            uniquemirnas.append(line['name'])
    uniquemirnas = set(uniquemirnas)
    
    # log the counts
    countsdf = np.log(countsdf + 1)
    
    manojset = "MP1 MP2 MP9 MP20 MP21 MP24 MP34 MP35 MP36 MP38 MP42.ACTG MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA".split()
    peterset1 = "PK11 PK21 PK24 PK31 PK41 PK42 PK51 PK52 PK54".split()
    peterset2 = "PK11 PK12 PK21 PK22 PK31 PK32 PK41 PK51 PK52 PK53".split()
    
    # print matrix
    countsdf.ix[uniquemirnas,manojset].to_csv("manoj_top50.csv", sep=",", header=True)
    countsdf.ix[uniquemirnas,peterset1].to_csv("peter1_top50.csv", sep=",", header=True)
    countsdf.ix[uniquemirnas,peterset2].to_csv("peter2_top50.csv", sep=",", header=True)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--mirnalist', nargs="+", help='miRNA list')
    p.add_argument('--counts', nargs="+", help='sample counts files')
    p.add_argument('--fullext', help="specify everything that comes after the sample name in the file name of the counts files")
    args = p.parse_args()
    
    main(args)