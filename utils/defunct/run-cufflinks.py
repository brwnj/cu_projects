#!/usr/bin/env python
# encoding: utf-8
"""
Runs Cufflinks on supplied mapped reads.

Created by Joe Brown on 2011-12-12.
"""
import argparse
import os
import sys
from subprocess import Popen


def runcufflinks(bam, geneannotation, referencegenome):
    bam = os.path.abspath(bam)
    annotation = os.path.abspath(geneannotation)
    genome = os.path.abspath(referencegenome)
    script = '%s/cufflinks.sh' % os.path.dirname(bam)
    sh = open(script, 'w')
    sh.write('#!/bin/sh\n')
    sh.write('#PBS -l walltime=72:00:00,nodes=1:ppn=8,mem=30gb\n')
    sh.write('#PBS -j oe\n')
    sh.write('cufflinks'\
                ' --output-dir %s'\
                ' --num-threads 8'\
                ' --GTF %s'\
                ' --multi-read-correct'\
                #' --library-type fr-secondstrand'\
                ' --upper-quartile-norm'\
                ' --label %s'\
                ' --frag-bias-correct %s'\
                ' %s' % (os.path.dirname(bam),
                         annotation,
                         bam.split("/")[-2],
                         genome, bam))
    sh.close()
    Popen(['qsub', script])


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group('required arguments')
    required.add_argument("-b", "--bam", dest="bam", nargs='+', 
                            help="Aligned reads (.bam). Using *.bam will also work.")
    required.add_argument("-g", "--gene-annotation", dest="gene", 
                            help="Gene annotation file (.gtf)")
    required.add_argument("-r", "--ref-sequence", dest="seq", 
                            help="Reference sequence (.fa)")
    args = p.parse_args()
    if args.bam and args.gene and args.seq:
        for b in args.bam:
            run_cufflinks(b, args.gene, args.seq)
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()