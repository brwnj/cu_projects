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

reps = "MP1 MP2 MP9 MP20 MP21 MP24 MP34 MP35 MP36 MP38 MP42.ACTG MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA".split()

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
        
        for rep in reps:
            casedata[rep][b.name] = {}
            tracklist = []
            intensities = chromosome[b.start:b.stop, ["%s.%s" % (rep, strand)]]
            casedata[rep][b.name] = get_peak_max(intensities)

caselist = casedata.keys()
features = casedata[caselist[0]].keys()

print "\t" + "\t".join(k for k in caselist)
for feature in features:
    print feature + "\t" + "\t".join(str(casedata[k][feature]) for k in caselist)