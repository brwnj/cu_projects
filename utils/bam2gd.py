#! /usr/bin/env python
# encoding: utf-8
"""
Converts bam to genomedata archive, saving stranded bedgraphs and bigwigs
along the way.
"""
import os
import sys
import os.path as op
from bsub import bsub
from itertools import izip

def bam2bg(bams, csizes):
    """given bams and chromosome sizes, writes two stranded bedgraph files."""
    sub = bsub("bam2bg", q="short", verbose=True)
    jobs = []
    for bam in bams:
        for symbol, strand in izip("- +".split(), "neg pos".split()):
            cmd = "genomeCoverageBed -bg -strand %s -ibam %s -g %s > %s_%s.bedgraph" \
                    % (symbol, bam, csizes, op.splitext(bam)[0], strand)
            jobs.append(sub(cmd))
    return jobs

def bg2bw(bgs, csizes):
    """convert bedgraphs to bigwigs using bedGraphToBigWig."""
    sub = bsub("bg2bw", q="short", verbose=True)
    for bg in bgs:
        sub("bedGraphToBigWig %s %s %s.bw" % (bg, csizes, op.splitext(bg)[0]))

def load_genome(bgs, fasta):
    """loads data into a genome archive with sequence data."""
    tracks = "".join(["-t %s=%s " \
                    % (op.splitext(op.basename(bg))[0], bg) for bg in bgs])
    cmd = "genomedata-load -v --directory-mode -s %s %s genomedata" \
                    % (fasta, tracks)
    print cmd
    return bsub("genomedata", verbose=True)(cmd)

def cleanup(bgs):
    """gzips the bedgraphs."""
    sub = bsub("zip", q="short", verbose=True)
    for bg in bgs:
        sub("gzip -f %s" % bg)

def main(args):
    # bam -> bedgraph
    bsub.poll(bam2bg(args.BAM, args.CHROM_SIZES))
    # bedgraph -> bw
    bgs = ["%s_%s.bedgraph" \
                    % (op.splitext(bam)[0], strand) \
                    for bam in args.BAM for strand in "neg pos".split()]
    bg2bw(bgs, args.CHROM_SIZES)
    # load genomedata
    bsub.poll(load_genome(bgs, args.FASTA))
    # gzip the bedgraphs    
    cleanup(bgs)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('CHROM_SIZES')
    p.add_argument('FASTA', help='reference for whole genome')
    p.add_argument('BAM', nargs="+",
                    help='bam(s) to convert to genomedata archive')
    args = p.parse_args()
    main(args)