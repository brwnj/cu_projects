#!/usr/bin/env python
# encoding: utf-8
"""
methods to aid in bulding a nextgen pipeline 
"""

import sys
import os
import os.path as op
import fnmatch
import pandas as pd
from bsub import bsub
from pybedtools import BedTool

def getfilelist(path, pattern):
    files = []
    for root, dirnames, filenames in os.walk(path):
          for filename in fnmatch.filter(filenames, pattern):
              files.append(op.join(root, filename))
    return files

def clobber_previous(results_path):
    import shutil
    shutil.rmtree(results_path)
    os.makedirs(results_path)

def fastqc(fastqc_path, samples, data_path, read_ext="fastq.gz"):
    sub = bsub("qc", verbose=True)
    for sample in samples:
        fastq = getfilelist(data_path, "%s.%s" % (sample, read_ext))
        assert(len(fastq) == 1)
        fastq = fastq[0]
        qcresult = "%s/%s_fastqc.zip" % (data_path.rstrip("/"), sample)
        if op.exists(qcresult): continue
        cmd = "%s --outdir %s %s" % (script, data_path, fastq)
        sub(cmd)

def counts(samples, result_path, peak_ext="_peaks.bed", bam_ext="bam"):
    # get the consensus peaks
    f = open("%s/peak_coordinates.bed" % result_path, 'w')
    x = BedTool()
    sub = bsub("counts", R="select[mem>16] rusage[mem=16] span[hosts=1]", 
                    verbose=True)
    consensus = x.multi_intersect(i=getfilelist(result_path, "*%s" % peak_ext))
    for c in consensus:
        replicate_counts = c.name
        if replicate_counts < 2: continue
        
        fields = [c.chrom, c.start, c.stop, "%s:%d-%d\n" % \
                    (c.chrom, c.start, c.stop)]
        f.write("\t".join(map(str, fields)))
    f.close()
    # get counts for each sample
    jobs = []
    countfiles = []
    for sample in samples:
        bams = getfilelist(result_path, "%s.%s" % (sample, bam_ext))
        assert(len(bams) == 1)
        outdir = result_path + "/" + sample
        countsresult = outdir + "/" + sample + ".counts"
        countfiles.append(countsresult)
        if op.exists(countsresult): continue
        cmd = "bedtools coverage -abam %s -b %s > %s" % \
                    (bams[0], f.name, countsresult)
        jobid = sub(cmd)
        jobs.append(jobid)
    bsub.poll(jobs)
    # counts to matrix
    allcounts = {}
    for cf in countfiles:
        cfname = op.basename(cf).split(".bam")[0]
        casecounts = {}
        for toks in reader(cf, header="chrom start stop name a_overlaps_in_b \
                    b_with_nonzero length_b frac_b_nonzero".split()):
            casecounts[toks['name']] = int(toks['a_overlaps_in_b'])
        allcounts[cfname] = casecounts
    countsdf = pd.DataFrame(allcounts)
    countsdf.to_csv(result_path + "/sample_counts.csv", sep=",", header=True)

def cleanup(path, extensions=['bed', 'xls']):
    for ext in extensions:
        for f in getfilelist(path, "*.%s" % ext):
            cmd = "gzip -f %s" % f
            bsub("zip", q="idle")(cmd)

def macs(samples, result_path, cmdstr):
    jobs = []
    for sample in samples:
        bams = getfilelist(result_path, "%s.bam" % sample)
        assert(len(bams) == 1)
        outdir = "%s/%s" % (result_path, sample)
        macsresult = "%s/%s_peaks.xls" % (outdir, sample)
        if op.exists(macsresult) or op.exists(macsresult + ".gz"): continue
        cmd = cmdstr.format(bams[0], sample)
        jobid = bsub("macs", cwd=outdir,
                        R="select[mem>16] rusage[mem=16] span[hosts=1]",
                        verbose=True)(cmd)
        jobs.append(jobid)
    return jobs

def alignment_stats(results_path, picard_path, ref_fasta):
    for bam in getfilelist(results_path, "*.bam"):
        cmd = "samtools index %s" % bam
        if not op.exists("%s.bai" % bam):
            jobid = bsub("index", verbose=True)(cmd)
            bsub.poll(jobid)
        cmd = "java -Xmx8g -jar %s/CollectMultipleMetrics.jar \
                INPUT=%s REFERENCE_SEQUENCE=%s ASSUME_SORTED=true OUTPUT=metrics \
                PROGRAM=CollectAlignmentSummaryMetrics \
                PROGRAM=QualityScoreDistribution \
                PROGRAM=MeanQualityByCycle" % (picard_path, bam, ref_fasta)
        bsub("alignment_summary", verbose=True)(cmd)

def gsnap(samples, reads_path, results_path, gmap_db, cmd_str, read_ext="bam"):
    """align reads for each sample according to the command string."""
    jobs = []
    sub = bsub("align", n="5", R="select[mem>28] rusage[mem=28] span[hosts=1]",
                    verbose=True)
    for sample in samples:
        fastqs = getfilelist(reads_path, "%s.%s" % (sample, read_ext))
        assert(len(fastqs) == 1)
        fastq = fastqs[0]
        
        out = "%s/%s" % (results_path, sample)
        if not op.exists(out):
            os.makedirs(out)
        
        align_result = "%s/%s.bam" % (out, sample)
        if op.exists(align_result): continue
        
        cmd = cmd_str.format(gmap_db, fastq, sample, align_result)
        jobs.append(sub(cmd))
    return jobs

def trim(path, pattern, trim_ext="trim.fastq.gz"):
    """uses https://github.com/lh3/seqtk"""
    jobs = []
    sub = bsub("trim", verbose=True)
    for fastq in getfilelist(path, pattern):
        trimresult = "%s.%s" % (fastq.split(".fastq", 1)[0], trim_ext)
        if op.exists(trimresult): continue
        cmd = "seqtk trimfq %s | gzip -c > %s" % (fastq, trimresult)
        jobs.append(sub(cmd))
    return jobs

def trim_adapter(datadir, adapters, in_file_ext="fastq.gz", out_file_ext="trim.fastq.gz"):
    """trim adapters using ea-utils"""
    jobs = []
    for fastq in getfilelist(datadir, "*.%s" % file_ext):
        trimresult = fastq.split(in_file_ext, 1)[0] + out_file_ext
        if op.exists(trimresult): continue
        cmd = "fastq-mcf " + adapters + " " + fastq + " | gzip -c > " + trimresult
        jobid = bsub("trim", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs