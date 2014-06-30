#!/usr/bin/env python
# encoding: utf-8
"""
Cambier
"""
from bsub import bsub
import os
import sys
import fnmatch

BASE = "/vol1/home/brownj/projects/cambier"
DATA = BASE + "/data"
RESULTS = BASE + "/results/common"
SAMPLES = "11 12 13 14 1 2 3 4 5 6 7 8".split()
INDEX="/vol1/home/brownj/ref/rum/mm9"

def getfilelist(datadir, searchpattern):
    """returns list of fastqs matching the search pattern, eg. 1_*."""
    files = []
    for root, dirnames, filenames in os.walk(datadir):
          for filename in fnmatch.filter(filenames, searchpattern):
              files.append(os.path.join(root, filename))
    return files


def getfilename(fname):
    """returns the file basename minus the extension."""
    return os.path.basename(fname).split(".", 1)[0]


def fastqc():
    """qc for single or paired-end data."""
    fastqc="/vol1/home/brownj/opt/fastqc/fastqc"
    for sample in SAMPLES:
        outdir = RESULTS + "/" + sample
        
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        
        for fastq in getfilelist(DATA, sample + "_*"):
            # see if result exists -- fastqc naming convention uses a portion of the read file name
            qcresult = getfilename(fastq) + "_fastqc.zip"
            if os.path.exists(outdir + "/" + qcresult): continue
            
            cmd = fastqc + " --outdir " + outdir + " --threads 4 " + fastq
            bsub("fastqc", verbose=True)(cmd)


def concat():
    """join reads from all lanes."""
    concatjobs = []
    for sample in SAMPLES:
        fastqs = getfilelist(DATA, sample + "_*")
        
        # check for output
        concatresult = DATA + "/" + sample + ".fq.gz"
        if os.path.exists(concatresult): continue
        
        assert(len(fastqs) == 2)
        
        cmd = "zcat " + " ".join(fastqs) + " | gzip -c > " + concatresult
        concatjobs.append(bsub("concat_reads", verbose=True)(cmd))


def align():
    """align reads using rum"""
    alignjobs = []
    for sample in SAMPLES:
        fastqs = getfilelist(sample + ".fq.gz")
        assert(len(fastqs) == 1)
        outdir = RESULTS + "/" + sample
        alignresult = outdir + "/" + sample + ".bam"
        if os.path.exists(alignresult): continue
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        cmd = "rum_runner align -v -i " + INDEX + " -o " + outdir + " --chunks 5 --name " + sample + " " + DATA + "/" + fastq[0]
        alignjobs.append(bsub("align_reads", "-n 5", verbose=True)(cmd))


def cleanup():
    """take care of the mess left by rum."""
    sortjobs = []
    for sam in getfilelist(RESULTS, "RUM.sam")
        outdir = os.path.dirname(sam)
        sample = outdir.rsplit("/", 1)[1]
        cmd = "samtools view -ShuF 4 " + sam + " | samtools sort -o - " + sample + ".temp -m 9500000000 > " + sample + ".bam"
        sortjobs.append(bsub("sam2bam", "-cwd " + outdir, verbose=True)(cmd))
        cmd = "gzip *.fa RUM_Unique RUM_NU"
        bsub("compress", "-cwd " + outdir, verbose=True)(cmd)
    return sortjobs


def indexbams():
    for sample in SAMPLES:
        for bam in getfilelist(RESULTS, sample + ".bam"):
            assert(len(bam) == 1)
            outdir = os.path.dirname(bam)
            cmd = "samtools index " + bam
            bsub("indexing", "-cwd " + outdir, verbose=True)(cmd)


def cuffdiff():
    """docstring for cuffdiff"""
    pass


def main():
	fastqc()
	# does this actually work?
	bsub.poll(concat())
	bsub.poll(align())
	bsub.poll(cleanup())
	indexbams()
    # cuffdiff()
	


if __name__ == '__main__':
	main()