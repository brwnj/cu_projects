#!/usr/bin/env python
# encoding: utf-8
"""
Walter
"""
from bsub import bsub
from subprocess import Popen, PIPE
import os
import os.path as op
import shlex
import sys
import fnmatch
import shutil

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


PI = "walter"
BASE = "/vol1/home/brownj/projects/walter"
DATA = BASE + "/data/20121005"
RESULTS = BASE + "/results/common/"
SAMPLES = "E1T1_Inf E1T1_Uninf E1T24_Inf E1T24_Uninf E1T2_Inf E1T2_Uninf \
            E1T8_Inf E1T8_Uninf E2T1_Inf E2T1_Uninf E2T24_Inf E2T24_Uninf \
            E2T2_Inf E2T2_Uninf E2T8_Inf E2T8_Uninf E3T1_Inf E3T1_Uninf \
            E3T24_Inf E3T24_Uninf E3T2_Inf E3T2_Uninf E3T8_Inf E3T8_Uninf".split()
INDEX = "/vol1/home/brownj/ref/rum/hg19"


def getfilelist(datadir, searchpattern):
    """assumes CASAVA naming convention (sample_index_etc_etc.fastq)"""
    files = []
    for root, dirnames, filenames in os.walk(datadir):
          for filename in fnmatch.filter(filenames, searchpattern):
              files.append(op.join(root, filename))
    return files


def getfilename(fname):
    """returns the basename of the reads file."""
    return op.basename(fname).split(".", 1)[0]


def extract(fname, out):
    """extracts the file in its current directory"""
    cmd = "zcat " + fname + " > " + out
    jobid = bsub("unzip", verbose=True)(cmd)
    return jobid


def fastqc(samples):
    """qc for single or paired-end data"""
    fastqc="/vol1/home/brownj/opt/fastqc/fastqc"
    for sample in samples:
        outdir = RESULTS + sample
        
        if not op.exists(outdir):
            os.makedirs(outdir)
        
        for fastq in getfilelist(DATA, sample + "_*"):
            # see if result exists -- fastqc naming convention uses a portion of the read file name
            qcresult = getfilename(fastq) + "_fastqc.zip"
            if op.exists(outdir + "/" + qcresult): continue
            
            cmd = fastqc + " --outdir " + outdir + " --threads 4 " + fastq
            bsub("fastqc", n="4", verbose=True)(cmd)


def trimadapters(datadir):
    """trim adapters using ea-utils"""
    jobs = []
    adapters = "/vol1/home/brownj/projects/walter/data/20121005/adapters.fa"
    for fastq in getfilelist(datadir, "*.fastq.gz"):
        trimresult = op.dirname(fastq) + "/" + op.basename(fastq).split(".fastq", 1)[0] + ".trm.fq.gz"
        if op.exists(trimresult): continue
        cmd = "fastq-mcf " + adapters + " " + fastq + " | gzip -c > " + trimresult
        jobid = bsub(PI + ".trimadapter", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def rumalign(samples, index, genome):
    """align to index using rum"""
    jobs = []
    for sample in samples:
        fastqs = getfilelist(DATA, sample + "_*.trm.fq.gz")
        assert(len(fastqs) == 1)
        
        outdir = RESULTS + sample
        alignresult = outdir + "/" + sample + "." + genome + ".bam"
        alternatealignresult = outdir + "/RUM.sam"
        if op.exists(alignresult) or op.exists(alternatealignresult): continue
        if not op.exists(outdir):
            os.makedirs(outdir)
        
        gzipfastq = fastqs[0]
        fastq = outdir + "/" + op.splitext(op.basename(gzipfastq))[0]
        if not op.exists(fastq):
            bsub.poll(extract(gzipfastq, fastq))
        #--limit-nu n
        # Limits the number of ambiguous mappers in the final output by removing all reads that map to n locations or more.
        cmd = "rum_runner align -v -i " + index + " -o " + outdir + " --chunks 5 --dna --name " + sample + " " + fastq
        jobid = bsub(PI + ".align_reads", n="5", R="select[mem>28] rusage[mem=28] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def cleanup(genome):
    """take care of the mess left by rum."""
    jobs = []
    for sam in getfilelist(RESULTS, "RUM.sam"):
        outdir = op.dirname(sam)
        try:
            [os.remove(fastq) for fastq in getfilelist(outdir, "*.fastq")]
        except OSError:
            # no fastq found
            pass
        sample = outdir.rsplit("/", 1)[1]
        
        cmd = "gzip *.fa RUM_Unique RUM_NU RUM_NU.cov RUM_Unique.cov"
        bsub("compress", cwd=outdir, verbose=True)(cmd)
        
        bam = outdir + "/" + sample + "." + genome + ".bam"
        if op.exists(bam): continue
        cmd = "samtools view -ShuF 4 " + sam + " | samtools sort -o - " + sample + ".temp -m 9500000000 > " + bam
        jobs.append(bsub(PI + ".distill_rum", cwd=outdir, verbose=True)(cmd))
    return jobs


def removesams(species):
    for sam in getfilelist(RESULTS, "RUM.sam"):
        outdir = op.dirname(sam)
        sample = outdir.rsplit("/", 1)[1]
        if op.exists(outdir + "/" + sample + "." + species + ".bam"):
            os.remove(sam)


def removefastqs():
    try:
        [os.remove(fastq) for fastq in getfilelist(RESULTS, "*.fastq")]
        [os.remove(fastq) for fastq in getfilelist(RESULTS, "*.fq")]
    except OSError:
        pass


def bowtiealign(samples, index, genome):
    """align to index using bowtie"""
    jobs = []
    for sample in samples:
        fastqs = getfilelist(DATA, sample + "_*.trm.fq.gz")
        # single end
        assert(len(fastqs) == 1)
        
        outdir = RESULTS + sample
        alignresult = outdir + "/" + sample + "." + genome + ".bam"
        if op.exists(alignresult):continue
        if not op.exists(outdir):
            os.makedirs(outdir)
        
        gzipfastq = fastqs[0]
        fastq = outdir + "/" + os.path.splitext(op.basename(gzipfastq))[0]
        if not op.exists(fastq):
            bsub.poll(extract(gzipfastq, fastq))
        
        cmd = "bowtie -p4 --best --sam -q " + index + " " + fastq + " | samtools view -ShuF4 - | samtools sort -o - " + sample + ".temp -m 9500000000 > " + alignresult
        jobid = bsub(PI + ".bowtie", n="4", R="select[mem>20] rusage[mem=20] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def macs(samples, resultsdir):
    """macs14 -t <bam> -f BAM -n <sample name> -g hs -w --single-profile"""
    jobs = []
    for sample in samples:
        bams = getfilelist(resultsdir, sample + "*.hg19.bam")
        assert(len(bams) == 1)
        outdir = resultsdir.rstrip("/") + "/" + sample
        macsresult = outdir + "/" + sample + "_peaks.xls"
        if op.exists(macsresult): continue
        cmd = "macs14 -t " + bams[0] + " -f BAM -n " + sample + " -g hs -w --single-profile"
        # writes to directory in which it was executed
        # output files include
            # E3T1_Uninf_MACS_wiggle/
            # E1T1_Inf_model.r
            # E1T1_Inf_peaks.bed
            # E1T1_Inf_peaks.xls
            # E1T1_Inf_summits.bed
            # E1T1_Uninf_model.r
            # E1T1_Uninf_peaks.bed
            # E1T1_Uninf_peaks.xls
            # E1T1_Uninf_summits.bed
        jobid = bsub(PI + ".macs", cwd=outdir, R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


# not working. giving up for now.

# def findconsensus():
#     """run in the working dir that ran macs
#     bedtools multiinter -cluster -i *.bed | awk 'BEGIN{OFS=FS="\t"}{print $1,$2,$3,"macs_peak_"NR}' > consensus.bed"""
#     bedfiles = getfilelist(os.getcwd(), "*_peaks.bed")
#     cmd = "bedtools multiinter -cluster -i " + " ".join([b for b in bedfiles])
#     print cmd
#     mergedbed = Popen(shlex.split(cmd), stdout=PIPE, stderr=PIPE, shell=False)
#     fout = open("consensus.bed", "wb")
#     for i, line in enumerate(mergedbed.stdout, start=1):
#         print line
#         # line = line.rstrip("\r\n").split("\t")
#         # chrom = line[0]
#         # start = line[1]
#         # stop = line[2]
#         # fields = [chrom, start, stop, "macs_peak_%d" % i]
#         # fout.write("\t".join(map(str, fields)))


def counts(samples, resultsdir):
    """get counts over peaks regions for each sample"""
    jobs = []
    # the merged peaks file
    consensus = getfilelist(BASE + "/results", "consensus.bed*")
    assert(len(consensus) == 1)
    consensus = consensus[0]
    for sample in samples:
        bams = getfilelist(resultsdir, sample + "*.hg19.bam")
        assert(len(bams) == 1)
        outdir = resultsdir.rstrip("/") + "/" + sample
        countsresult = outdir + "/" + sample + ".counts"
        if op.exists(countsresult): continue
        cmd = "bedtools coverage -abam " + bams[0] + " -b " + consensus + " > " + countsresult
        jobid = bsub(PI + ".counts", R="select[mem>16] rusage[mem=16] span[hosts=1]", verbose=True)(cmd)
        jobs.append(jobid)
    return jobs


def main():
    hairpinindex = "/vol1/home/brownj/ref/mirbase/19/hairpin19"
    matureindex = "/vol1/home/brownj/ref/mirbase/19/mature19"
    tuberculosisindex = "/vol1/home/brownj/ref/tuberculosis/H37Rv"
    fastqc(SAMPLES)
    bsub.poll(trimadapters(DATA))
    
    # Rum
    bsub.poll(rumalign(SAMPLES, INDEX, "hg19"))
    bsub.poll(cleanup("hg19"))
    removesams("hg19")
            
    # Bowtie
    bsub.poll(bowtiealign(SAMPLES, matureindex, "mature"))
    bsub.poll(bowtiealign(SAMPLES, tuberculosisindex, "H37Rv"))
    removefastqs()
    
    # MACS (peak calling)
    bsub.poll(macs(SAMPLES, RESULTS))
    # use bedtools merge to create a consensus of peaks
    # i did this manually using the command in the docstring of the following method
    # bsub.poll(findconsensus())
    bsub.poll(counts(SAMPLES, RESULTS))


if __name__ == '__main__':
    main()