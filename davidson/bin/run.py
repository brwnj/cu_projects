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
import sys
import fnmatch
import shutil


def getfilelist(datadir, searchpattern):
    """assumes CASAVA naming convention (sample_index_etc_etc.fastq)"""
    files = []
    for root, dirnames, filenames in os.walk(datadir):
          for filename in fnmatch.filter(filenames, searchpattern):
              files.append(os.path.join(root, filename))
    return files


def extract(fname, out):
    """extracts the file in its current directory"""
    cmd = "zcat " + fname + " > " + out
    jobid = bsub("unzip", verbose=True)(cmd)
    return jobid


def fastqc(samples, resultsdir, datadir):
    """qc for single or paired-end data"""
    fastqc="/vol1/home/brownj/opt/fastqc/fastqc"
    for sample in samples:
        outdir = resultsdir + sample
        
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        for fastq in getfilelist(datadir, sample + "_*.fastq.gz"):
            # gives full path to data directory, need results/basename
            qcresult = os.path.splitext(fastq)[0].rstrip(".fastq") + "_fastqc.zip"
            
            if os.path.exists(outdir + "/" + os.path.basename(qcresult)): continue
            
            cmd = fastqc + " --outdir " + outdir + " --threads 4 " + fastq
            bsub("fastqc", verbose=True)(cmd)


def trim(samples, datadir):
    """writes to DATA. uses seqtk."""
    jobs = []
    for sample in samples:
        for fastq in getfilelist(datadir, sample + "_*.fastq.gz"):
            trimresult = fastq.split(".fastq", 1)[0] + ".trm.fq.gz"
            if os.path.exists(trimresult): continue
            
            # https://github.com/lh3/seqtk
            # trim low qual bases from both ends
            cmd = "seqtk trimfq " + fastq + " | gzip -c > " + trimresult
            jobid = bsub("seqtk", verbose=True)(cmd)
            jobs.append(jobid)
    return jobs


def join(samples, datadir, basedir):
    """joins paired-end data into SSAKE format."""
    jobs = []
    for sample in samples:
        # sort for ordering: R1 then R2
        fastqs = sorted(getfilelist(datadir, sample + "_*.trm.fq.gz"))
        
        # check for output
        joinresult = datadir + "/" + sample + ".jnd.fa.gz"
        if os.path.exists(joinresult) or os.path.exists(joinresult + ".gz"): continue
        
        assert(len(fastqs) == 2)
        
        # usage: join_reads.py R1 R2
        cmd = "python " + basedir + "/bin/join_reads.py " + " ".join(fastqs) + " | gzip -c > " + joinresult
        jobid = bsub("join_reads", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def assemble(samples, datadir, resultsdir):
    """assemble using iSSAKE."""
    jobs = []
    for sample in samples:
        fastas = getfilelist(datadir, sample + ".jnd.fa.gz")
        assert(len(fastas) == 1)
        
        gzipfasta = fastas[0]
        
        outdir = resultsdir + sample
        
        #contigresult = getfilelist(resultsdir, sample + "_*.contigs")
        #if len(contigresult) > 1: continue
        
        fasta = outdir + "/" + os.path.splitext(os.path.basename(gzipfasta))[0]
        if not os.path.exists(fasta):
            bsub.poll(extract(gzipfasta, fasta))
        
        cmd = "iSSAKE -f " + fasta + " -m 15 -o 2 -r 0.7 -t 0 -c -b " + sample + " -z 50 -p 1 -v 1 -d 200 -e 0.75 -k 2 -a 0.7"
        jobid = bsub("assemble_reads", R="select[mem>20] rusage[mem=20] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def removefastas(resultsdir):
    try:
        [os.remove(fasta) for fasta in getfilelist(resultsdir, "*.fasta")]
    except OSError:
        pass


def main():
    base = "/vol1/home/brownj/projects/davidson"
    data = base + "/data/20120924"
    results = base + "/results/common/"
    samples = ["1","2","3","4","5","6"]
    
    fastqc(samples, base + "/results/common/", data)
    bsub.poll(trim(samples, data))
    bsub.poll(join(samples, data, base))
    bsub.poll(assemble(samples, data, results))
    #atexit.register
    removefastas(base + "/results/common")


if __name__ == '__main__':
    main()