#!/usr/bin/env python
# encoding: utf-8
"""
Given fasta, gff, and bam, parses for TA, annotates feature, and find coverage.
"""
import re
import sys
import tempfile
from cogent.parse.fasta import MinimalFastaParser
from toolshed import nopen
from subprocess import Popen, PIPE
from progressbar import ProgressBar, Percentage, Bar

__author__ = "Joe Brown"
__author_email__ = "brwnjm@gmail.com"


def getlocation(gff, chrom, start, stop):
    # overlap this TA to see if it's in a gene
    # assumes each TA only overlaps one gene
    # if overlap exists, location is updated
    ta = open(tempfile.mktemp(suffix=".bed"), 'w')
    ta.write('%s\t%d\t%d\n' % (chrom, start, stop))
    ta.close()

    p = Popen(["bedtools", "intersect", "-wb", "-a", ta.name, "-b", gff], stdout=PIPE, shell=False)
    location = 0
    genename = "."
    for e in p.stdout:
        e = e.strip("\r\n").split("\t")
        genestart = int(e[6])
        genestop = int(e[7])
        geneattrs = e[11]
        genename = re.findall(r'Name=([\w\.]+)', geneattrs)[0]
        location = 100 * ((start - genestart) / (genestop - genestart - 1.))
    p.wait()
    return genename, "%.0f%%" % location


def main(args):
    # i already know the approx. number of TA sites
    pbar = ProgressBar(widgets=[Percentage(), Bar()], maxval=92000).start()
    counter = 0
    # every fasta entry
    for f in MinimalFastaParser(nopen(args.fasta)):
        chrom = f[0]
        sequence = f[1]
        match = [m.start() for m in re.finditer('TA', sequence)]
        
        # every match in individual fasta entry
        for m in match:
            name = "TA_%d" % counter
            start = m
            stop = m + 2
            coverage = 0
            
            genename, location = getlocation(args.gff, chrom, start, stop)
            
            # count the number of reads starting at the TA site
            samfile = Popen(["samtools", "view", args.bam, "%s:%d-%d" % (chrom, start, stop)], stdout=PIPE, shell=False)
            for sam in samfile.stdout:
                readstart = int(sam.split("\t")[3])
                if readstart == start + 1:
                    coverage += 1
            samfile.wait()
            
            # for naming the TAs
            counter =+ 1
            # update the progress bar
            pbar.update(counter)
            
            fields = (chrom, start, stop, name, genename, location, coverage)
            print "\t".join(map(str, fields))
    pbar.finish()


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument('fasta')
    p.add_argument('gff')
    p.add_argument('bam')
    args = p.parse_args()

    main(args)