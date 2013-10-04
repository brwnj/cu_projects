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
fisher_script = "/vol1/home/brownj/devel/polya/fisher_test.py"
fisher_results = "/vol1/home/brownj/projects/polya/results/common/fisher_results"
metadata = "/vol1/home/brownj/projects/polya/results/common/hub/metadata.tsv"
sample_column = "alias"
comparison_column = "comparisons"
files = glob("/vol1/home/brownj/projects/polya/results/common/*/*.counts.txt.gz")
fisher_submit = bsub("fisher", q="short", P="pillai_kabos_polya", verbose=True)

def comparison_complete(results_dir, result_file):
    for f in os.listdir(results_dir):
        if os.path.splitext(os.path.basename(result_file))[0] in os.path.basename(f):
            return True
    return False

for t in reader(metadata, header=True):
    for compareto, strand in zip(t[comparison_column].split(","), ['pos', 'neg']):
        if len(compareto) == 0: continue
        a = t[sample_column]
        b = compareto
        fisher_result = "{fisher_results}/{a}_to_{b}.{strand}.fisher.txt.gz".format(**locals())

        if comparison_complete(fisher_results, fisher_result):
            print >>sys.stderr, ">> fisher comparison complete for", a, "and", b
            continue
        print >>sys.stderr, ">> fisher: comparing", a, "and", b

        file_a = "{results}/{a}/{a}.{strand}.counts.txt.gz".format(**locals())
        file_b = "{results}/{b}/{b}.{strand}.counts.txt.gz".format(**locals())
        assert os.path.exists(file_a) and os.path.exists(file_b)

        fisher_cmd = ("python {fisher_script} {file_a} {file_b} "
                        "| gzip -c > {fisher_result}").format(**locals())
        fisher_submit(fisher_cmd)
