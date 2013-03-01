#! /usr/bin/env python
# encoding: utf-8
"""
Converts bam(s) to genomedata archive, saving stranded bedgraphs and bigwigs
along the way. File names are parsed to create the stranded track names.
"""
import os
import sys
import fnmatch
import os.path as op
from bsub import bsub
from itertools import izip

def extract(fname, out):
    """unzips fname to out."""
    sub = bsub("unzip", q=short, verbose=True)
    cmd = "zcat " + fname + " > " + out
    return sub(cmd)

def bam2bg(bams, csizes, clobber):
    """given bams and chromosome sizes, writes two stranded bedgraph files."""
    sub = bsub("bam2bg", q="short", verbose=True)
    jobs = []
    for bam in bams:
        for symbol, strand in izip("- +".split(), "neg pos".split()):
            bg = "%s_%s.bedgraph" % (op.splitext(bam)[0], strand)
            if op.exists(bg) and not clobber: continue
            if op.exists("%s.gz" % bg) and not clobber:
                jobs.append(extract(bg))
                continue
            cmd = "genomeCoverageBed -bg -strand %s -ibam %s -g %s > %s" \
                    % (symbol, bam, csizes, bg)
            jobs.append(sub(cmd))
    return jobs

def bg2bw(bgs, csizes, clobber):
    """convert bedgraphs to bigwigs using bedGraphToBigWig."""
    sub = bsub("bg2bw", q="short", verbose=True)
    for bg in bgs:
        bw = "%s.bw" % op.splitext(bg)[0]
        if op.exists(bw) and not clobber: continue
        sub("bedGraphToBigWig %s %s %s" % (bg, csizes, bw))

def load_genome(bgs, fasta_or_dir, out):
    """loads data into a genome archive with sequence data."""
    sub = bsub("genomedata", verbose=True)
    # reference files
    fastas = []
    if op.isdir(fasta_or_dir):
        for root, dirnames, filenames in os.walk(fasta_or_dir):
              # allows for .fasta, .fa, and .gz variants
              for filename in fnmatch.filter(filenames, "*.fa*"):
                  fastas.append(op.join(root, filename))
    else:
        fastas.append(fasta_or_dir)
    # track names
    tracks = "".join(["-t %s=%s " \
                    % (op.splitext(op.basename(bg))[0], bg) for bg in bgs])
    # job ids in order to wait until completion
    jobs = []
    for fasta in fastas:
        cmd = "genomedata-load -v --directory-mode -s %s %s %s" \
                        % (fasta, tracks, out)
        jobs.append(sub(cmd))
    return jobs

def cleanup(bgs):
    """gzips the bedgraphs."""
    sub = bsub("zip", q="short", verbose=True)
    for bg in bgs:
        sub("gzip -f %s" % bg)

def main(args):
    # bam -> bedgraph
    bsub.poll(bam2bg(args.BAM, args.CHROM_SIZES, args.clobber))
    # bedgraph -> bw
    bgs = ["%s_%s.bedgraph" \
                    % (op.splitext(bam)[0], strand) \
                    for bam in args.BAM for strand in "neg pos".split()]
    bg2bw(bgs, args.CHROM_SIZES, args.clobber)
    # bedgraphs -> genomedata
    bsub.poll(load_genome(bgs, args.FASTA, args.output))
    # gzip the bedgraphs
    cleanup(bgs)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
                    formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('CHROM_SIZES', help="chromosome sizes for genome")
    p.add_argument('FASTA', help='fasta or folder containing fasta(s)')
    p.add_argument('BAM', nargs="+",
                    help='bam(s) to convert to genomedata archive')
    p.add_argument('-o', '--output', default="genomedata",
                    help="genomedata archive name [ %(default)s ]")
    p.add_argument('--clobber', action="store_true",
                    help="overwrite existing intermediate files [ %(default)s ]")
    main(p.parse_args())