#!/usr/bin/env python
# encoding: utf-8
"""
Tophat.py
When dealing with paired-end data, set the flag and pass the pairs, 1 then 2, 
into the script. It will account for more than one run at a time, and will
skip any trailing reads without a pair.

Created by Joe Brown on 2011-12-16.
"""
import argparse
import os
import sys
from subprocess import Popen


def run_tophat(files, index, annotation="", pairedend=False, mateinnerdistance="", matestddev=""):
    """Change the Tophat parameters inline.
    files = .fastq
    index = genome.fa or hg19all.fa
    annotation = gene annotation file (.gtf)
    pairedend = library type
    mateinnerdistance = fragment length minus (read length * 2)
    matestddev = fragment length std dev
    """
    sideone = ""
    for i,f in enumerate(files):
        label = os.path.basename(f).split("_")[0].lower()
        f = os.path.abspath(f)
        read_dir = os.path.dirname(f)
        output_dir = '%s/%s' % (read_dir, label)
        if pairedend:
            if i % 2 is 0:
                sideone = "%s" % f
                continue
            else:
                f = "%s %s" % (sideone, f)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        index = os.path.abspath(index).rstrip(".fa")
        script = '%s/tophat.sh' % output_dir
        sh = open(script, 'w')
        sh.write("#!/bin/sh\n")
        sh.write("#PBS -l walltime=240:00:00,nodes=1:ppn=8,mem=30gb\n")
        sh.write("#PBS -j oe\n")
        sh.write("cd %s\n" % output_dir)
        if pairedend:
            sh.write("tophat\
            --num-threads 8 \
            --segment-mismatches 2 \
            --mate-inner-dist %s \
            --mate-std-dev %s \
            --segment-length 21 \
            --no-novel-juncs \
            --transcriptome-index /mnt/storage3/brownj/reference/Homo_sapiens/Ensembl/GRCh37/transcriptome_data/known \
            --transcriptome-only \
            --output-dir %s \
            %s %s\n" % (mateinnerdistance, matestddev, output_dir, index, f))
        else:
            sh.write("tophat --num-threads 8 --segment-mismatches 3 --output-dir %s %s %s\n" % (output_dir, index, f))
        sh.write("mv accepted_hits.bam RNASEQ_%s_.bam" % label)
        #sh.write("gzip %s" % f) #f.split(" ") for paired-end
        sh.close()
        Popen(['qsub', script])
    sys.exit(0)


def main():
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group("required arguments")
    required.add_argument("-i","--index", dest="index",
            help="Genome reference index.")
    required.add_argument("-f","--files", dest="reads", nargs='+', 
            help="FASTQ files.")
    pairedend = p.add_argument_group('paired-end arguments')
    pairedend.add_argument("-p","--paired-end", dest="pairedend", action='store_true', default=False, 
            help="Set for paired-end data. See description for details.")
    pairedend.add_argument("--mate-inner-distance",dest='mateinnerdistance', type=str,
            help="The expected inner distance between pairs per Tophat specifications.")
    pairedend.add_argument("--mate-std-dev", dest='matestddev', type=str,
            help="The standard deviation for distribution of inner distances between pairs.")
    pairedend.add_argument("-a","--gene-annotation", dest='geneannotation', default="", 
            help='Gene annotation file.')
    # hard-coded at the moment
    #pairedend.add_argument('-t','--transcriptome-index',dest='transcriptome',help='Transcriptome index.')
    args = p.parse_args()
    if args.reads and args.index:
        if args.pairedend:
            run_tophat(args.reads, args.index, args.geneannotation, \
                        args.pairedend, args.mateinnerdistance, \
                        args.matestddev)
        else:
            run_tophat(args.reads, args.index)
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()