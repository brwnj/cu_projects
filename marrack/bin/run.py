#!/usr/bin/env python
# encoding: utf-8
"""
Marrack
"""
from bsub import bsub
import os
import os.path as op
import sys
import fnmatch
import shutil


def getfilelist(datadir, searchpattern):
    files = []
    for root, dirnames, filenames in os.walk(datadir):
          for filename in fnmatch.filter(filenames, searchpattern):
              files.append(op.join(root, filename))
    return files


def extract(fname, out):
    """extracts the file in its current directory"""
    cmd = "zcat " + fname + " > " + out
    jobid = bsub("unzip", verbose=True)(cmd)
    return jobid


def fastqc(samples, datadir, resultsdir):
    """qc for single or paired-end data"""
    fastqc="/vol1/home/brownj/opt/fastqc/fastqc"
    for sample in samples:
        fastqs = getfilelist(datadir, sample + ".fastq.gz")
        assert(len(fastqs) == 1)
        
        outdir = resultsdir.rstrip("/") + "/" + sample
        
        if not op.exists(outdir):
            os.makedirs(outdir)
        
        
        qcresult = outdir + "/" + sample + "_fastqc.zip"
        if op.exists(qcresult): continue
        cmd = fastqc + " --outdir " + outdir + " " + fastqs[0]
        bsub("fastqc", verbose=True)(cmd)


def trim(path, pattern):
    """uses https://github.com/lh3/seqtk"""
    jobs = []
    for fastq in getfilelist(path, pattern):
        trimresult = fastq.split(".fastq", 1)[0] + ".trim.fastq.gz"
        if os.path.exists(trimresult): continue

        cmd = "seqtk trimfq " + fastq + " | gzip -c > " + trimresult
        jobs.append(bsub("seqtk", verbose=True)(cmd))
    return jobs


def rum(samples, datadir, resultsdir, index):
    """align to index using rum"""
    jobs = []
    for sample in samples:
        fastqs = getfilelist(datadir, sample + ".trim.fastq.gz")
        assert(len(fastqs) == 1)
        
        outdir = resultsdir + "/" + sample
        alignresult = outdir + "/" + sample + ".bam"
        alternatealignresult = outdir + "/RUM.sam"
        if op.exists(alignresult) or op.exists(alternatealignresult): continue
        
        gzipfastq = fastqs[0]
        fastq = outdir + "/" + op.splitext(op.basename(gzipfastq))[0]
        if not op.exists(fastq):
            bsub.poll(extract(gzipfastq, fastq))

        cmd = "rum_runner align -v -i " + index + " -o " + outdir + " --chunks 5 --dna --nu-limit 2 --variable-length-reads --name " + sample + " " + fastq
        jobid = bsub("rum", n="5", R="select[mem>28] rusage[mem=28] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def postprocessrum(resultsdir):
    """take care of the mess left by rum."""
    jobs = []
    for sam in getfilelist(resultsdir, "RUM.sam"):
        outdir = op.dirname(sam)
        try:
            [os.remove(fastq) for fastq in getfilelist(outdir, "*.fastq")]
        except OSError:
            pass
        sample = outdir.rsplit("/", 1)[1]
        
        cmd = "gzip -f *.fa RUM_Unique RUM_NU RUM_NU.cov RUM_Unique.cov"
        bsub("postprocessrum", q="idle", cwd=outdir, verbose=True)(cmd)
        
        bam = outdir + "/" + sample + ".bam"
        if op.exists(bam): continue
        cmd = "samtools view -ShuF 4 " + sam + " | samtools sort -o - " + sample + ".temp -m 9500000000 > " + bam
        jobs.append(bsub("postprocessrum", cwd=outdir, verbose=True)(cmd))
    return jobs


def macs(samples, resultsdir, control):
    jobs = []
    for sample in samples:
        bams = getfilelist(resultsdir, sample + ".bam")
        assert(len(bams) == 1)
        
        # control
        if bams[0].contains(control): continue
        controlbam = getfilelist(resultsdir, control + ".bam")
        assert(len(controlbam) == 1)
        
        outdir = resultsdir.rstrip("/") + "/" + sample
        macsresult = outdir + "/" + sample + "_peaks.xls"
        
        if op.exists(macsresult) or op.exists(macsresult + ".gz"): continue
        
        cmd = "macs14 -t " + bams[0] + "-c " + controlbam[0] + " -f BAM -n " + sample + " -g mm -w --single-profile"
        jobid = bsub("macs", cwd=outdir, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def cleanup(path):
    """it'd be a good idea to not run this on the data dir"""
    exts = ['bed', 'xls']
    for ext in exts:
        for f in getfilelist(path, "*." + ext):
            cmd = "gzip -f " + f
            bsub("zip", q="idle")(cmd)
    try:
        [os.remove(sam) for sam in getfilelist(path, "*.sam")]
    except OSError:
        pass


def main():
    samples = ['RS_input_CCGTCC_L005_R1_001',
                'RS_iso_ATGTCA_L005_R1_001',
                'RS_tbet_CTTGTA_L005_R1_001']
    control = 'RS_input_CCGTCC_L005_R1_001'
    datadir = "/vol1/home/brownj/projects/marrack/data/20121101"
    resultsdir = "/vol1/home/brownj/projects/marrack/results/common"
    rumindex = "/vol1/home/brownj/ref/rum/mm9"
    
    fastqc(samples, datadir, resultsdir)
    bsub.poll(trim(datadir, "*R1_001.fastq.gz"))
    bsub.poll(rum(samples, datadir, resultsdir, rumindex))
    bsub.poll(postprocessrum(resultsdir))
    bsub.poll(macs(samples, resultsdir, control))
    cleanup(resultsdir)


if __name__ == '__main__':
    main()