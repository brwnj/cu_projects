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


def macs(samples, resultsdir, controls):
    jobs = []
    for sample in samples:
        bams = getfilelist(resultsdir, sample + ".bam")
        assert(len(bams) == 1)
        
        # control
        if "Input" in bams[0]: continue
        controlbam = ""
        for control in controls:
            if control.split("_")[0] == bams[0].split("_")[0]:
                controlbam = getfilelist(resultsdir, control + ".bam")[0]
        
        outdir = resultsdir + "/" + sample
        macsresult = outdir + "/" + sample + "_peaks.xls"
        
        if op.exists(macsresult) or op.exists(macsresult + ".gz"): continue
        
        cmd = "macs14 -t " + bams[0] + " -c " + controlbam + " -f BAM -n " + sample + " -g hs -w --single-profile"
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
        bams = getfilelist(resultsdir, sample + "*.bam")
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
        cfname = op.basename(cf).split(".bam")[0]
        casecounts = {}
        for toks in reader(cf, header="chrom start stop name a_overlaps_in_b b_with_nonzero length_b frac_b_nonzero".split()):
            casecounts[toks['name']] = int(toks['a_overlaps_in_b'])
        allcounts[cfname] = casecounts
    countsdf = pd.DataFrame(allcounts)
    countsdf.to_csv(resultsdir + "/sample_counts.csv", sep=",", header=True)


def main():
    samples = ['2Som_chip1_GCCAAT_L006_R1_001',
               '2Som_chip2_GTCCGC_L006_R1_001',
               '2Som_Input_GTGAAA_L006_R1_001',
               '31hpt_Chip1_CAGATC_L006_R1_001',
               '31hpt_Chip2_ACAGTG_L006_R1_001',
               '31hpt_Input_TGACCA_L006_R1_001']
    controls = ['2Som_Input_GTGAAA_L006_R1_001',
                '31hpt_Input_TGACCA_L006_R1_001']
    datadir = "/vol1/home/brownj/projects/artinger/data/20121101"
    resultsdir = "/vol1/home/brownj/projects/artinger/results/common"
    rumindex = "/vol1/home/brownj/ref/rum/hg19"
    
    fastqc(samples, datadir, resultsdir)
    bsub.poll(trim(datadir, "*R1_001.fastq.gz"))
    bsub.poll(rum(samples, datadir, resultsdir, rumindex))
    bsub.poll(postprocessrum(resultsdir))
    bsub.poll(macs(samples, resultsdir, controls))
    cleanup(resultsdir)
    counts(samples, resultsdir)


if __name__ == '__main__':
    main()