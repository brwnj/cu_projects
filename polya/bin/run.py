#!/usr/bin/env python
# encoding: utf-8
"""
MP51 control
MP52 treat1
MP53 treat2
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


def extract(fin, fout):
    """extract fin to fout"""
    cmd = "zcat " + fin + " > " + fout
    jobid = bsub("unzip", verbose=True)(cmd)
    return jobid


def scarftofastq(scarf, fastq):
    """
    Converts scarf to gzipped fastq.
    Returns job id.
    """
    script = "/vol1/home/brownj/projects/polya/bin/scarftofastq.py"
    cmd = "python " + script + " " + scarf + " | gzip -c > " + fastq
    return bsub("scarftofastq", verbose=True)(cmd)


def removefastas(resultsdir):
    try:
        [os.remove(fasta) for fasta in getfilelist(resultsdir, "*.fasta")]
    except OSError:
        pass


def bowtie(fastq, index):
    """Align reads using bowtie."""
    sample = os.path.splitext(os.path.basename(fastq))[0]
    cmd = "bowtie -p4 --best --sam -q " + index + " " + fastq + " | samtools view -ShuF4 - | samtools sort -o - " + sample + ".temp -m 9500000000 > " + sample + ".bam"
    return bsub("bowtie", n="4", R="select[mem>20] rusage[mem=20] span[hosts=1]", verbose=True)(cmd)


def main():
    fqdir = "/vol1/home/brownj/projects/polya/data/20121008/"
    outdir = "vol1/home/brownj/projects/polya/results/common/"
    samples = ["MPG1"]
    
    #bsub.poll(scarftofastq("/vol1/home/brownj/projects/polya/data/20121008/s_1_sequence.txt", "/vol1/home/brownj/projects/polya/data/20121008/MPG1.fastq.gz"))
    #bsub.poll(demultiplex("/vol1/home/brownj/projects/polya/data/20121008/MPG1.fastq.gz", "/vol1/home/brownj/projects/polya/data/20121008/MPG1.dmplx.fastq.gz"))
    # bowtie("/vol1/home/brownj/projects/polya/data/20121008/MPG1.fastq", "/vol3/home/jhessel/projects/bowtie/indices/hg18")
    bowtie("/vol1/home/brownj/projects/polya/data/20121008/MPG1.ACTG.fastq", "/vol3/home/jhessel/projects/bowtie/indices/hg18")


if __name__ == '__main__':
    main()