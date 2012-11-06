#!/usr/bin/env python
# encoding: utf-8
"""
Table of snoRNA to case showing the max read intensity.
"""
from pybedtools import BedTool
from genomedata import Genome
from numpy import nansum
from collections import defaultdict
import sys

def get_peak_max(intensities):
    """Finds the maximum peak intensity value from ndarray."""
    maximum = 0
    for intensity in intensities:
        # sum intensity across samples
        sum_intensity = nansum(intensity)
        if sum_intensity > maximum:
            maximum = sum_intensity
    return maximum


cases = {'BT474herc':"PK23",
         'MCF7':"PK11 PK31 PK51",
         'MCF7estr':"PK12 PK32",
         'MDA231':"PK24 PK42 PK54",
         'BT474':"PK21 PK41 PK52",
         'BT474estr':"PK22 PK53",
         'HS27A':"MP1 MP21 MP35",
         'HS5':"MP2 MP20 MP34",
         'hMSC':"MP36 MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA",
         'BMEC':"MP42.ACTG MP45.ACTG MP45.TCGA",
         'HUVEC':"MP24 MP38"}


snorna = BedTool("/vol1/home/brownj/projects/hits-clip/data/20120815/snoRNA.bed.gz")

casedata = defaultdict(dict)
with Genome("/vol1/home/brownj/projects/hits-clip/data/combined.genomedata") as genome:    
    
    for b in snorna:
        chromosome = genome[b.chrom]
        
        strand = ""
        if b.strand is "+":
            strand = "pos"
        else:
            strand = "neg"
        
        for case, replicates in cases.iteritems():
            casedata[case][b.name] = {}
            tracklist = []
            for r in replicates.split():
                tracklist.append("%s.%s" % (r, strand))
            intensities = chromosome[b.start:b.stop, tracklist]            
            casedata[case][b.name] = get_peak_max(intensities)

caselist = casedata.keys()
features = casedata[caselist[0]].keys()

print "#snoRNA\t" + "\t".join(k for k in caselist)    
for feature in features:
    print feature + "\t" + "\t".join(str(casedata[c][feature]) for c in caselist)