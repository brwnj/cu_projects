#!/usr/bin/env python
# encoding: utf-8
"""
Duval
"""
from bsub import bsub
from pybedtools import BedTool
from toolshed import reader
import os
import os.path as op
import pandas as pd
import sys
import fnmatch
import shutil

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def getfilelist(datadir, searchpattern):
    """assumes CASAVA naming convention (sample_index_etc_etc.fastq)"""
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


def solid2fastq(samples, datadir):
    script = "/vol2/home/brentp/src/bfast-git/scripts/solid2fastq"
    jobs = []
    for sample in samples:
        csfastas = getfilelist(datadir, sample + "*.csfasta.gz")
        quals = getfilelist(datadir, sample + "*.qual.gz")
        assert(len(csfastas) == 1)
        assert(len(quals) == 1)
        if op.exists(datadir + "/" + sample + ".fastq.gz"): continue
        cmd = script + " -z -Z -o " + sample + " " + csfastas[0] + " " + quals[0]
        jobid = bsub("solid2fastq", cwd=datadir, verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


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


def bowtiealign(samples, datadir, resultsdir, index, genome):
    """align to index using bowtie"""
    jobs = []
    for sample in samples:
        fasta = datadir + "/" + sample + ".csfasta"
        qual = datadir + "/" + sample + ".qual"
        
        outdir = resultsdir.rstrip("/") + "/" + sample
        alignresult = outdir + "/" + sample + "." + genome + ".bam"
        if op.exists(alignresult): continue
        
        cmd = "bowtie -p4 -m1 -v1 -f -C --best --strata --chunkmbs 512 --trim3 25 --sam " + index + " -Q " + qual + " " + fasta + " | samtools view -ShuF4 - | samtools sort -o - " + sample + ".temp -m 9500000000 > " + alignresult        
        jobid = bsub("bowtie", n="4", R="select[mem>20] rusage[mem=20] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def novoalign(samples, datadir, resultsdir, index, genome):
    jobs = []
    for sample in samples:
        fastqs = getfilelist(datadir, sample + ".fastq.gz")
        assert(len(fastqs) == 1)
        
        outdir = resultsdir.rstrip("/") + "/" + sample
        alignresult = outdir + "/" + sample + "." + genome + ".bam"
        if op.exists(alignresult): continue
        if not op.exists(outdir):
            os.makedirs(outdir)
        gzipfastq = fastqs[0]
        fastq = outdir + "/" + op.splitext(op.basename(gzipfastq))[0]
        if not op.exists(fastq):
            bsub.poll(extract(gzipfastq, fastq))
        cmd = "novoalignCS -c 1 -d " + index + " -f " + fastq + " -F BFASTQ -o SAM -r Random -e 100 -s 8 -l 20 | samtools view -ShuF4 - | samtools sort -o - " + sample + ".temp -m 9500000000 > " + alignresult
        jobid = bsub("novoalign", n="1", R="select[mem>20] rusage[mem=20] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def macs(samples, resultsdir):
    jobs = []
    for sample in samples:
        bams = getfilelist(resultsdir, sample + ".hg19.bam")
        assert(len(bams) == 1)
        outdir = resultsdir.rstrip("/") + "/" + sample
        macsresult = outdir + "/" + sample + "_peaks.xls"
        if op.exists(macsresult) or op.exists(macsresult + ".gz"): continue
        cmd = "macs14 -t " + bams[0] + " -f BAM -n " + sample + " -g hs -w --single-profile"
        # writes to directory in which it was executed
        jobid = bsub("macs", cwd=outdir, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def cleanup(path):
    """it'd be a good idea not to run this on the data dir"""
    exts = ['bed', 'xls']
    for ext in exts:
        for f in getfilelist(path, "*." + ext):
            cmd = "gzip -f " + f
            bsub("zip", q="idle")(cmd)
    try:
        [os.remove(fastq) for fastq in getfilelist(path, "*.fastq")]
        [os.remove(fastq) for fastq in getfilelist(path, "*.fq")]
        [os.remove(fastq) for fastq in getfilelist(path, "*.csfasta")]
        [os.remove(fastq) for fastq in getfilelist(path, "*.qual")]
    except OSError:
        pass


def counts(samples, resultsdir):
    """docstring"""
    # get the consensus peaks
    f = open(resultsdir + "/peak_coordinates.bed", 'w')
    x = BedTool()
    consensus = x.multi_intersect(i=getfilelist(resultsdir, "*peaks.bed.gz"))
    for c in consensus:
        replicate_counts = c.name
        if replicate_counts < 2: continue
        
        fields = [c.chrom, c.start, c.stop, "%s:%d-%d\n" % (c.chrom, c.start, c.stop)]
        f.write("\t".join(map(str, fields)))
    f.close()
    
    # get counts for each sample
    jobs = []
    countfiles = []
    for sample in samples:
        bams = getfilelist(resultsdir, sample + "*.hg19_novoalign.bam")
        assert(len(bams) == 1)
        outdir = resultsdir.rstrip("/") + "/" + sample
        countsresult = outdir + "/" + sample + ".counts"
        countfiles.append(countsresult)
        if op.exists(countsresult): continue
        cmd = "bedtools coverage -abam " + bams[0] + " -b " + f.name + " > " + countsresult
        jobid = bsub(sample + "_counts", R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    bsub.poll(jobs)
    
    # counts to matrix
    allcounts = {}
    for cf in countfiles:
        cfname = op.basename(cf).split(".hg19_novoalign.bam")[0]
        casecounts = {}
        for toks in reader(cf, header="chrom start stop name a_overlaps_in_b b_with_nonzero length_b frac_b_nonzero".split()):
            casecounts[toks['name']] = int(toks['a_overlaps_in_b'])
        allcounts[cfname] = casecounts
    countsdf = pd.DataFrame(allcounts)
    countsdf.to_csv(resultsdir + "/sample_counts.csv", sep=",", header=True)


def trimreads(samples, datadir):
    script = "/vol1/home/brownj/projects/duval/bin/solid-trimmer.py"
    jobs = []
    for sample in samples:
        fastas = getfilelist(datadir, sample + ".csfasta.gz")
        quals = getfilelist(datadir, sample + ".qual.gz")
        # single end
        assert(len(fastas) == 1)
        assert(len(quals) == 1)
        
        trimresult = datadir + "/" + sample + ".trm_F3.csfasta"
        prefix = sample + ".trm"
        if op.exists(trimresult) or op.exists(trimresult + ".gz"): continue

        fasta = fastas[0]
        qual = quals[0]
        
        cmd = "python " + script + " -c " + fasta + " -q " + qual + " -p " + prefix + " --max-ns 3 --moving-average 7:12 --min-read-length 18"
        jobid = bsub("trimming", cwd=datadir, R="select[mem>8] rusage[mem=8] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def main():
    samples = ["Hela_Pit1_T220A_0", "Hela_Pit1_T220A_1", "Hela_Pit1_T220A_2",
                "Hela_Pit1_T220D_0", "Hela_Pit1_T220D_1", "Hela_Pit1_T220D_2",
                "Hela_Pit1_WT_0", "Hela_Pit1_WT_1", "Hela_Pit1_WT_2",
                "Hela_Pit_T220E_0", "Hela_Pit_T220E_1", "Hela_Pit_T220E_2",
                "HelaWT_IgG", "HelaWT_NoAb"]
    btindex = "/vol1/home/brownj/ref/hg19/hg19_c"
    novoindex = "/vol1/home/brownj/ref/hg19/hg19"
    datadir = "/vol1/home/brownj/projects/duval/data/20121015"
    resultsdir = "/vol1/home/brownj/projects/duval/results/common"
    
    # bsub.poll(solid2fastq(samples, datadir))
    # fastqc(samples, datadir, resultsdir)
    # bsub.poll(trimreads(samples, datadir))
    # bsub.poll(bowtiealign(samples, datadir, resultsdir, btindex, "hg19"))
    # bsub.poll(novoalign(samples, datadir, resultsdir, novoindex, "hg19_novoalign"))
    # bsub.poll(macs(samples, resultsdir))
    # cleanup(resultsdir)
    # consensus peaks
    # bedtools coverage
    # data -> matrix
    counts(samples, resultsdir)


if __name__ == '__main__':
    main()