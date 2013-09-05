#!/usr/bin/env python
# encoding: utf-8
"""
Only the comparisons defined in the metadata sheet will be performed.
"""
import os
import sys
import tempfile
import subprocess
from bsub import bsub
from glob import glob
from random import randint
from toolshed import reader

results = "/vol1/home/brownj/projects/polya/results/common"
dexseq_script = "/vol1/home/brownj/devel/polya/run_dexseq.R"
dexseq_results = "/vol1/home/brownj/projects/polya/results/common/dexseq_results"
metadata = "/vol1/home/brownj/projects/polya/results/common/hub/metadata.tsv"
sample_column = "alias"
comparison_column = "comparisons"
files = glob("/vol1/home/brownj/projects/polya/results/common/*/*.counts.txt.gz")
submit = bsub("dexseq", P="pillai_kabos_polya", q="normal", n="4", R="span[hosts=1]", verbose=True)

def replicate(fname, n=5):
    tmp = open(tempfile.mkstemp(suffix=".txt", dir=".")[1], 'w')
    for c in reader(fname, header=['name','count']):
        count = int(c['count'])
        low = count - n if count - n > 0 else 0
        count = randint(low, count + n)
        tmp.write("{name}\t{count}\n".format(name=c['name'], count=count))
    tmp.close()
    return tmp.name

for t in reader(metadata, header=True):
    for compareto, strand in zip(t[comparison_column].split(","), ['pos', 'neg']):
        if len(compareto) == 0: continue
        comparison_complete = False
        a = t[sample_column]
        b = compareto
        result = "{dexseq_results}/{a}_vs_{b}.{strand}.dexseq.txt".format(**locals())
        for f in os.listdir(dexseq_results):
            if os.path.splitext(os.path.basename(result))[0] in os.path.basename(f):
                comparison_complete = True
        if comparison_complete:
            print >>sys.stderr, ">> comparison complete for", a, "and", b
            continue
        print >>sys.stderr, ">> comparing", a, "and", b
        file_a = "{results}/{a}/{a}.{strand}.counts.txt.gz".format(**locals())
        file_b = "{results}/{b}/{b}.{strand}.counts.txt.gz".format(**locals())
        assert os.path.exists(file_a) and os.path.exists(file_b)
        rep_a = replicate(file_a)
        rep_b = replicate(file_b)
        cmd = ("Rscript {script} {a},{a}x {file_a},{rep_a} "
                "{b},{b}x {file_b},{rep_b} "
                "{result}").format(script=dexseq_script, a=a, file_a=file_a,
                                    rep_a=rep_a, b=b, file_b=file_b,
                                    rep_b=rep_b, result=result)
        submit(cmd)
