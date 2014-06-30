#!/usr/bin/env python
# encoding: utf-8
"""
Artinger
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

def trim(path, pattern):
    """uses https://github.com/lh3/seqtk"""
    jobs = []
    for fastq in getfilelist(path, pattern):
        trimresult = "%s.trim.fastq.gz" % fastq.split(".fastq", 1)[0]
        if op.exists(trimresult): continue
        cmd = "seqtk trimfq %s | gzip -c > %s" % (fastq, trimresult)
        jobs.append(bsub("trim", verbose=True)(cmd))
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
        
        cmd = cmd_str.format(gmap_db, fastq, sample, align_result)
        jobid = bsub("align", n="5", R="select[mem>28] rusage[mem=28] span[hosts=1]", verbose=True)(cmd)
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

def macs(samples, resultsdir, controls, cmdstr):
    jobs = []
    for sample in samples:
        bams = getfilelist(resultsdir, sample + ".bam")
        assert(len(bams) == 1)

        # control
        if "Input" in bams[0]: continue
        # controlbam = ""
        # for control in controls:
        #     sampleid = control.split("_")[0]
        #     if sampleid in bams[0]:
        #         controlbam = getfilelist(resultsdir, control + ".bam")[0]
        
        outdir = resultsdir + "/" + sample
        macsresult = outdir + "/" + sample + "_peaks.xls"
        
        if op.exists(macsresult) or op.exists(macsresult + ".gz"): continue
        
        cmd = cmdstr.format(bams[0], sample)
        print cmd
        jobid = bsub("macs", cwd=outdir, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs

def cleanup(path):
    exts = ['bed', 'xls']
    for ext in exts:
        for f in getfilelist(path, "*." + ext):
            cmd = "gzip -f " + f
            bsub("zip", q="idle")(cmd)

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
    fastqc_script="/vol1/home/brownj/opt/fastqc/fastqc"
    picard = "/vol1/home/brownj/opt/picard-tools-1.79"
    rumindex = "/vol1/home/brownj/ref/rum/zebrafish"
    reference_fasta = "/vol1/home/brownj/ref/zebrafish/Danio_rerio.Zv9.68.fa"
    gmapdb = "/vol1/home/brownj/ref/gmapdb"
    
    rumcmd = "rum_runner align -v -i %s -o {} --chunks 5 --dna --nu-limit 2 --variable-length-reads --name {} {}" % rumindex
    macscmd = "macs14 -t {} -f BAM -n {} -g 1400000000 -w --single-profile --call-subpeaks"
    gsnapcmd = "gsnap -D {} -d zebrafish --gunzip --npaths=1 --quiet-if-excessive --batch=5 --nofails --nthreads=4 --format=sam {} | samtools view -ShuF 4 - | samtools sort -o - {}.temp -m 9500000000 > {}"
    
    if args.clobber:
        clobber_previous(resultsdir)
    
    fastqc(fastqc_script, samples, datadir)
    bsub.poll(trim(datadir, "*R1_001.fastq.gz"))
    bsub.poll(gsnap(samples, datadir, resultsdir, gmapdb, gsnapcmd))
    # alignment_stats(resultsdir, picard, reference_fasta)
    bsub.poll(macs(samples, resultsdir, controls, macscmd))
    cleanup(resultsdir)
    # counts(samples, resultsdir)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clobber", action="store_true", help="clear all previous results")
    args = p.parse_args()
    main(args)