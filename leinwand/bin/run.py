#!/usr/bin/env python
# encoding: utf-8
"""
Leinwand
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

def getfilelist(path, pattern):
    files = []
    for root, dirnames, filenames in os.walk(path):
          for filename in fnmatch.filter(filenames, pattern):
              files.append(op.join(root, filename))
    return files

def extract(fin, fout):
    """extract fin to fout"""
    cmd = "zcat %s > %s" % (fin, fout)
    jobid = bsub("unzip", verbose=True)(cmd)
    return jobid

def fastqc(script, samples, data_path):
    for sample in samples:
        fastq = getfilelist(data_path, sample + ".fastq.gz")
        assert(len(fastq) == 1)
        fastq = fastq[0]
        qcresult = "%s/%s_fastqc.zip" % (data_path.rstrip("/"), sample)
        if op.exists(qcresult): continue
        cmd = "%s --outdir %s %s" % (script, data_path, fastq)
        bsub("qc", verbose=True)(cmd)

def trimadapter(datadir, adapters):
    """trim adapters using ea-utils"""
    jobs = []
    for fastq in getfilelist(datadir, "*.fastq.gz"):
        trimresult = op.dirname(fastq) + "/" + op.basename(fastq).split(".fastq", 1)[0] + ".trim.fastq.gz"
        if op.exists(trimresult): continue
        cmd = "fastq-mcf " + adapters + " " + fastq + " | gzip -c > " + trimresult
        jobid = bsub("trim", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs

def gsnap(samples, reads_path, results_path, gmap_db, cmd_str):
    """align reads for each sample according to the command string."""
    jobs = []
    for sample in samples:
        fastqs = getfilelist(reads_path, sample + ".trim.fastq.gz")
        assert(len(fastqs) == 1)
        fastq = fastqs[0]
        
        out = "%s/%s" % (results_path, sample)
        if not op.exists(out):
            os.makedirs(out)
        
        align_result = "%s/%s.bam" % (out, sample)
        if op.exists(align_result): continue
        
        cmd = cmd_str.format(gmap_db, fastq, sample, sample)
        jobid = bsub("align", n="5", R="select[mem>28] rusage[mem=28] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs

def alignment_stats(results_path, picard_path, ref_fasta):
    for bam in getfilelist(results_path, "*.bam"):
        stats_result = "%s_stats.txt" % op.splitext(bam)[0]
        if op.exists("%s.bai" % bam) and op.exists(stats_result): continue
        cmd = "samtools index %s" % bam
        jobid = bsub("index", verbose=True)(cmd)
        bsub.poll(jobid)
        cmd = "java -Xmx8g -jar %s/CollectMultipleMetrics.jar \
                INPUT=%s REFERENCE_SEQUENCE=%s ASSUME_SORTED=true OUTPUT=%s \
                PROGRAM=CollectAlignmentSummaryMetrics \
                PROGRAM=QualityScoreDistribution \
                PROGRAM=MeanQualityByCycle" % (picard_path, bam, ref_fasta, stats_result)
        bsub("alignment_summary", verbose=True)(cmd)

def macs(samples, resultsdir, cmdstr):
    jobs = []
    for sample in samples:
        bams = getfilelist(resultsdir, sample + ".bam")
        assert(len(bams) == 1)
        
        outdir = resultsdir + "/" + sample
        macsresult = outdir + "/" + sample + "_peaks.xls"
        
        if op.exists(macsresult) or op.exists(macsresult + ".gz"): continue
        
        cmd = cmdstr.format(bams[0], sample)
        jobid = bsub("macs", cwd=outdir, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs

def cleanup(path):
    exts = ['bed', 'xls']
    for ext in exts:
        for f in getfilelist(path, "*." + ext):
            cmd = "gzip -f " + f
            bsub("zip", q="idle")(cmd)
    try:
        [os.remove(sam) for sam in getfilelist(path, "*.sam")]
    except OSError:
        pass

def counts(samples, result_path):
    # get the consensus peaks
    f = open(result_path + "/peak_coordinates.bed", 'w')
    x = BedTool()
    consensus = x.multi_intersect(i=getfilelist(result_path, "*_peaks.bed"))
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
        bams = getfilelist(result_path, sample + "*.bam")
        assert(len(bams) == 1)
        outdir = result_path.rstrip("/") + "/" + sample
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
    countsdf.to_csv(result_path + "/sample_counts.csv", sep=",", header=True)

def clobber_previous(results_path):
    import shutil
    shutil.rmtree(results_path)
    os.makedirs(results_path)

def main(args):
    samples = ["MDX_22_AGTTCC_L003_R1_001",
               "MDX_23_ATGTCA_L003_R1_001",
               "MDX_24_CCGTCC_L003_R1_001",
               "WT_21_AGTCAA_L003_R1_001",
               "WT_25_GTAGAG_L003_R1_001",
               "WT_42_GTCCGC_L003_R1_001"]
    datadir = "/vol1/home/brownj/projects/leinwand/data/20121101"
    adapters = "%s/adapters.fa" % datadir
    resultsdir = "/vol1/home/brownj/projects/leinwand/results/common"
    fastqc_script="/vol1/home/brownj/opt/fastqc/fastqc"
    rumindex = "/vol1/home/brownj/ref/rum/mm9"
    
    rumcmd = "rum_runner align -v -i %s -o {} --chunks 5 --dna --nu-limit 2 --variable-length-reads --name {} {}" % rumindex
    macscmd = "macs14 -t {} -f BAM -n {} -g mm -w --single-profile"
    
    if args.clobber:
        clobber_previous(resultsdir)
    
    fastqc(fastqc_script, samples, datadir)
    bsub.poll(trimadapter(datadir, adapters))
    bsub.poll(rum(samples, datadir, resultsdir, rumcmd))
    bsub.poll(postprocessrum(resultsdir))
    bsub.poll(macs(samples, resultsdir, controls, macscmd))
    cleanup(resultsdir)
    counts(samples, resultsdir)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clobber", action="store_true", help="clear all previous results")
    args = p.parse_args()
    main(args)