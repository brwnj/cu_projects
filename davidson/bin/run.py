#!/usr/bin/env python
# encoding: utf-8
"""
Davidson

sample key:
1   z1+z5         a
2   z1+z5         b
3   z9+z19        a
4   z9+z19        b
5   b9+23ee+JK    a
6   b9+23ee+JK    b

"""
from bsub import bsub
import os
import os.path as op
import sys

sys.path.append('/vol1/home/brownj/projects/utils')
import ngseq

def join(samples, datadir, basedir):
    """joins paired-end data into SSAKE format."""
    jobs = []
    for sample in samples:
        # sort for ordering: R1 then R2
        fastqs = sorted(getfilelist(datadir, sample + "_*.trm.fq.gz"))
        
        # check for output
        joinresult = datadir + "/" + sample + ".jnd.fa.gz"
        if op.exists(joinresult) or op.exists(joinresult + ".gz"): continue
        
        assert(len(fastqs) == 2)
        
        # usage: join_reads.py R1 R2
        cmd = "python " + basedir + "/bin/join_reads.py " + " ".join(fastqs) + " | gzip -c > " + joinresult
        jobid = bsub("join_reads", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs

def assemble(samples, data_dir, results_dir, seed_fa):
    """assemble using SSAKE."""
    jobs = []
    sub = bsub("assemble_reads", R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)
    for sample in samples:
        fastas = getfilelist(datadir, sample + ".jnd.fa.gz")
        assert(len(fastas) == 1)
        
        gzipfasta = fastas[0]
        outdir = "%s/%s" % (results_dir, sample)
        
        fasta = outdir + "/" + op.splitext(op.basename(gzipfasta))[0]
        if not op.exists(fasta):
            bsub.poll(ngseq.extract(gzipfasta, fasta))
        
        cmd = "SSAKE -f " + fasta + " -s " + seed_fa + " -m 40 -o 50 -r 0.8 -b " + sample + " -p 1 -v 1 -d 200 -e 0.75 -k 10 -a 0.5 -x 50"
        jobs.append(sub(cmd))
    return jobs

def main():
    base = "/vol1/home/brownj/projects/davidson"
    data = base + "/data/20120924"
    results = base + "/results/common"
    samples = ["1","2","3","4","5","6"]
    beta = "path to trbv.fa"
    alpha = "path to trav.fa"
    
    # fastqc(samples, base + "/results/common/", data)
    # bsub.poll(trim(samples, data))
    # bsub.poll(join(samples, data, base))
    assemble(samples, data, results, beta)
    assemble(samples, data, results, alpha)
    # removefastas(base + "/results/common")

if __name__ == '__main__':
    main()