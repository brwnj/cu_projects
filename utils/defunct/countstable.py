#!/usr/bin/env python
# encoding: utf-8
import argparse
import os
import sys
from subprocess import Popen, PIPE
"""
Takes bams OR already processed bed files and creates a count matrix.
"""

def rawcounts(bam, referencebed):
    """Writes an annotated bed with read overlap information. The pertinent 
    information is column 12.
    This method will attempt to wait for jobs to complete from the compute
    nodes so it would not be a bad idea to run in the background by ending 
    your command with '&'.
    
    bam = aligned reads of which you want to measure overlaps
    referencebed = bed12 of genes created using getbed.py
    """
    allcountsfiles = []
    jobsrunning = {}
    for b in bam:
        b = os.path.abspath(b)
        referencebed = os.path.abspath(referencebed)
        folder = os.path.dirname(b)
        # assumes bam naming follows illumina read file naming convention
        label = os.path.basename(b).split("_")[1]
        countfile = '%s/%s.counts.bed' % (folder, label)
        allcountsfiles.append(countfile)
        script = '%s/%s.deseq.sh' % (folder, label)
        sh = open(script, 'w')
        sh.write('#!/bin/sh\n')
        sh.write('#PBS -j oe\n')
        # samtools view -bh -F 0x2 %s |
        sh.write("coverageBed -abam %s -b %s > %s"
                    % (b, referencebed, countfile))
        sh.close()
        if len(bam) > 1:
            jobsrunning[label] = Popen(['qsub', '-I', 'sh', script])
        else:
            Popen(['qsub', script])
    if len(bam) > 1:
        for k, v in jobsrunning.iteritems():
            jobsrunning[k].wait()
            sys.stderr.write("Job complete for sample: %s\n" % k)
        makecountstable(allcountsfiles)


def makecountstable(files):
    """Prints count data to stdout."""
    #TODO: needs to account for technical replicates by summing their counts into one column
    sys.stderr.write("Preparing counts table.\n")
    all_data = {}
    for f in files:
        # ie abspath/LABEL.counts.bed
        sample = os.path.basename(f).split(".")[0]
        all_data[sample] = {}
        for line in open(f):
            toks = line.split("\t")
            name = toks[3]
            #column verification
            try:
                count = str(int(toks[12]))
            except IndexError:
                sys.stderr.write("Incomplete count data in: %s\n" % f)
                sys.stderr.write("Error first encountered here:\n%s\n" % line)
                sys.exit(1)
            except ValueError:
                sys.stderr.write("Reference likely not bed12 format.\n")
                sys.stderr.write("Culprit file: %s\n" % f)
                sys.stderr.write("Error first encountered here:\n%s\n" % line)
                sys.exit(1)
            all_data[sample][name] = count
    samples = sorted(all_data)
    genes = all_data[samples[1]].keys()
    print "#ID\t" + "\t".join(samples)
    for gene in genes:
        print gene + "\t" + "\t".join(all_data[s][gene] for s in samples)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group("required arguments")
    required.add_argument("-b", "--bam", dest='bam', nargs='+', 
        help="Bam or bams of mapped reads. Space delimited list or *.bam.")
    required.add_argument("-a", "--annotation", dest='annotation', 
        help="Reference annotation file.")
    p.add_argument("--beds", dest="beds", nargs='+', 
        help="Annotated beds that have already been through coverageBed.")
    args = p.parse_args()
    if args.annotation and args.bam:
        rawcounts(args.bam, args.annotation)
    elif args.beds:
        makecountstable(args.beds)
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()