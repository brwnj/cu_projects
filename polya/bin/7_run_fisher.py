#!/usr/bin/env python
# encoding: utf-8
"""
Submit jobs for Fisher testing.
"""
import os
from bsub import bsub
from string import Template
from itertools import combinations

results = "/vol1/home/brownj/projects/polya/results/common"
projid = "pillai_kabos_polya"
script = "/vol1/home/brownj/devel/polya/fisher_test.py"

# get the samples
pksamples = [s for s in os.listdir(results) if s.startswith("PK")]
mpsamples = [s for s in os.listdir(results) if s.startswith("MP")]
# iterate through all combinations of MP and PK samples, submitting jobs
submit = bsub("fisher", q="short", P=projid, verbose=True)
out_template = Template("$results/fisher_results/${a}_to_$b.$strand.fisher.txt.gz")
cmd_template = Template("python $script $results/$a/$a.$strand.counts.txt.gz $results/$b/$b.$strand.counts.txt.gz | gzip -c > $out")
for strand in ["pos", "neg"]:
    # lazy
    for (a, b) in combinations(pksamples, 2):
        # write files to fisher_results
        out = out_template.substitute(**locals())
        cmd = cmd_template.substitute(**locals())
        submit(cmd)
    for (a, b) in combinations(mpsamples, 2):
        out = out_template.substitute(**locals())
        cmd = cmd_template.substitute(**locals())
        submit(cmd)
