#!/usr/bin/env python
# encoding: utf-8
"""
transfer miseq files to project/data/common
"""
import os
import re
import fnmatch
import pandas as pd
import subprocess as sp
from shutil import copyfile
from collections import defaultdict

def wget(prefix, url):
    sp.call("wget -r -l1 --no-parent -P {prefix} -A.fastq.gz {url}".format(prefix=prefix, url=url), shell=True)

data = "/vol1/home/brownj/projects/bennett/data/common"
meta_tsv = "/vol1/home/brownj/projects/bennett/data/common/metadata.tsv"
meta = pd.read_table(meta_tsv, index_col=['patient', 'barcode'])
run_from_path = re.compile(r"Data/(.*)/Data")
meta['run'] = meta.url.apply(lambda x: run_from_path.findall(x)[0])

# transfer from miseq to tesla
seen = set()
for (url, run) in zip(meta.url, meta.run):
    if run in seen: continue
    seen.add(run)
    run_folder = "{data}/genomics-cluster.ucdenver.pvt/MiSeq/Data/{run}".format(data=data, run=run)
    if os.path.exists(run_folder): continue
    # transfer these files from the MiSeq
    wget(data, url)

# list fastqs by run
fastqs = defaultdict(list)
for root, dirnames, filenames in os.walk("{data}/genomics-cluster.ucdenver.pvt".format(data=data)):
    for filename in fnmatch.filter(filenames, "*.fastq.gz"):
        path = os.path.join(root, filename)
        run = run_from_path.findall(path)[0]
        fastqs[run].append(path)

# copy over files and rename
for index, row in meta.iterrows():
    patient, barcode = index
    run_fastqs = fastqs[row['run']]
    index_read = "{data}/{patient}_{cell_type}_I1.fastq.gz".format(data=data, patient=patient, cell_type=row['cell_type'])
    r1 = "{data}/{patient}_{cell_type}_R1.fastq.gz".format(data=data, patient=patient, cell_type=row['cell_type'])
    r2 = "{data}/{patient}_{cell_type}_R2.fastq.gz".format(data=data, patient=patient, cell_type=row['cell_type'])

    for fastq in run_fastqs:
        base = os.path.basename(fastq)
        if not base.startswith(barcode): continue
        if not os.path.exists(index_read) and "I1" in base:
            copyfile(fastq, index_read)
        elif not os.path.exists(r1) and "R1" in base:
            copyfile(fastq, r1)
        elif not os.path.exists(r2) and "R2" in base:
            copyfile(fastq, r2)
        else:
            continue
