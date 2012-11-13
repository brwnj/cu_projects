#!/usr/bin/env python
# encoding: utf-8
"""
parallelize miso on lsf

specify run_miso.py script, miso index, bam, read length, overhang length, and output directory
"""
import sys
import fnmatch
import os
import os.path as op
from bsub import bsub


def getfilelist(path, pattern):
    files = []
    for root, dirnames, filenames in os.walk(path):
          for filename in fnmatch.filter(filenames, pattern):
              files.append(op.join(root, filename))
    return files


def main(args):
    if not op.exists(args.out):
        os.makedirs(args.out)
    
    jobids = []
    genelist = getfilelist(args.index, "*.pickle")
    for i, gene in enumerate(genelist):
        gname = op.splitext(op.basename(gene))[0]
        cmd = "python %s --read-len %s --overhang-len %s --settings-filename %s --compute-gene-psi %s %s %s %s" % (args.miso_script, args.read_length, args.overhang_length, args.miso_settings, gname, gene, args.bam, args.out)
        bsub("miso_" + gname, q=args.queue_name)(cmd, job_cap=args.job_cap)


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group('required')
    required.add_argument('--bam', '-b', required=True, help="reads file")
    required.add_argument('--index', '-x', required=True, help="miso index")
    required.add_argument('--out', '-o', required=True, help="sample output location")
    p.add_argument('--queue_name', default="normal", help="desired queue [ normal ]")
    p.add_argument('--job_cap', default=100, type=int, help="number of jobs to queue at a time [ 100 ]")
    p.add_argument('--overhang_length', default="1", help="overhang length [ 1 ]")
    p.add_argument('--read_length', default="100", help="read length [ 100 ]")
    p.add_argument('--miso_script', default="/vol1/software/modules-python/python/2.7.2/lib/python2.7/site-packages/misopy-0.4.4-py2.7-linux-x86_64.egg/misopy/run_miso.py", help="path to run_miso.py [ system version ]")
    p.add_argument('--miso_settings', default="/vol1/software/modules-python/python/2.7.2/lib/python2.7/site-packages/misopy-0.4.6-py2.7-linux-x86_64.egg/misopy/settings/miso_settings.txt", help="path to miso settings [ system version ]")
    args = p.parse_args()
    
    main(args)