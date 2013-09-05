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
from string import Template
from toolshed import reader

results = "/vol1/home/brownj/projects/polya/results/common"
dexseq_script = "/vol1/home/brownj/devel/polya/run_dexseq.R"
dexseq_results = "/vol1/home/brownj/projects/polya/results/common/dexseq_results"
fisher_script = "/vol1/home/brownj/devel/polya/fisher_test.py"
fishser_results = "/vol1/home/brownj/projects/polya/results/common/fisher_results"
metadata = "/vol1/home/brownj/projects/polya/results/common/hub/metadata.tsv"
sample_column = "alias"
comparison_column = "comparisons"
files = glob("/vol1/home/brownj/projects/polya/results/common/*/*.counts.txt.gz")
dexseq_submit = bsub("dexseq", P="pillai_kabos_polya", n="4", R="span[hosts=1]", verbose=True)
fisher_submit = bsub("fisher", q="short", P="pillai_kabos_polya", verbose=True)

def replicate(fname, n=5):
    """Generate DEXSeq replicate."""
    tmp = open(tempfile.mkstemp(suffix=".txt", dir=".")[1], 'w')
    for c in reader(fname, header=['name','count']):
        count = int(c['count'])
        low = count - n if count - n > 0 else 0
        count = randint(low, count + n)
        tmp.write("{name}\t{count}\n".format(name=c['name'], count=count))
    tmp.close()
    return tmp.name

def comparison_complete(results_dir, result_file):
    for f in os.listdir(results_dir):
        if os.path.splitext(os.path.basename(result_file))[0] in os.path.basename(f):
            return True
    return False

for t in reader(metadata, header=True):
    for compareto, strand in zip(t[comparison_column].split(","), ['pos', 'neg']):
        if len(compareto) == 0: continue
        comparison_complete = False
        a = t[sample_column]
        b = compareto
        dexseq_result = "{dexseq_results}/{a}_vs_{b}.{strand}.dexseq.txt".format(**locals())
        fisher_result = "{fisher_results}/{a}_vs_{b}.{strand}.fisher.txt.gz".format(**locals())

        if comparison_complete(dexseq_results, dexseq_result):
            print >>sys.stderr, ">> dexseq comparison complete for", a, "and", b
            continue
        print >>sys.stderr, ">> dexseq: comparing", a, "and", b
        
        if comparison_complete(fisher_results, fisher_result):
            print >>sys.stderr, ">> fisher comparison complete for", a, "and", b
            continue
        print >>sys.stderr, ">> fisher: comparing", a, "and", b
        
        file_a = "{results}/{a}/{a}.{strand}.counts.txt.gz".format(**locals())
        file_b = "{results}/{b}/{b}.{strand}.counts.txt.gz".format(**locals())
        assert os.path.exists(file_a) and os.path.exists(file_b)
        
        # create dexseq replicates
        rep_a = replicate(file_a)
        rep_b = replicate(file_b)
        dexseq_cmd = ("Rscript {script} {a},{a}x {file_a},{rep_a} "
                        "{b},{b}x {file_b},{rep_b} "
                        "{result}").format(script=dexseq_script, a=a,
                                    file_a=file_a, rep_a=rep_a, b=b,
                                    file_b=file_b, rep_b=rep_b, result=result)
        dexseq_submit(dexseq_cmd)
        
        fisher_cmd = ("python {fisher_script} {file_a} {file_b} | "
                        "gzip -c > {fisher_result}").format(**locals())
        fisher_submit(fisher_cmd)
