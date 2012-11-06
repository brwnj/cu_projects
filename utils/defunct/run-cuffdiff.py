#!/usr/bin/env python
# encoding: utf-8
"""
Runs Cuffdiff over samples specified in the Config file supplied to Cuffmerge.

Created by Joe Brown on 2011-12-20.
"""
import argparse
import os
import sys
import gzip
import bz2
import urllib
from itertools import izip


def nopen(f, mode="rb"):
    if not isinstance(f, basestring):
        return f
    if f.startswith("|"):
        p = Popen(f[1:], stdout=PIPE, stdin=PIPE, shell=True)
        if mode[0] == "r": return p.stdout
        return p
    return {"r": sys.stdin, "w": sys.stdout}[mode[0]] if f == "-" \
         else gzip.open(f, mode) if f.endswith((".gz", ".Z", ".z")) \
         else bz2.BZ2File(f, mode) if f.endswith((".bz", ".bz2", ".bzip2")) \
         else urllib.urlopen(f) if f.startswith(("http://", "https://",
             "ftp://")) \
        else open(f, mode)

    
def reader(fname, header=True, sep="\t"):
    line_gen = (l.rstrip("\r\n").split(sep) for l in nopen(fname))
    if header == True:
        header = line_gen.next()
        header[0] = header[0].lstrip("#")

    if header:
        for toks in line_gen:
            yield dict(izip(header, toks))
    else:
        for toks in line_gen:
            yield toks


def cuffdiff():
    """Reads the manifest file from Cuffmerge and runs samples through 
    in groups set by the manifest.
    """
    outputdir
    mergedgtf
    bowtieindex
    pathtobamsforgroup1
    pathtobamsforgroup2
    sh.write("#!/bin/sh\n")
    sh.write("#PBS -l walltime=48:00:00,nodes=1:ppn=8,mem=30gb\n")
    sh.write("#PBS -j oe\n")
    sh.write("cuffdiff \
                --output-dir $out \
                --labels AJ,LM \
                --num-threads 8 \
                --upper-quartile-norm \
                --total-hits-norm \
                --multi-read-correct \
                --library-type fr-secondstrand \
                --frag-bias-correct $ref \
                $mergedgtf? \
                group1 \
                group2")


def get_args(config):
    """Parses the config file and gets relevant arguments."""
    d_args = defaultdict(list)
    for l in reader(config, header=False):
        if 'index' in l[0]:
            d_args['index'] = os.path.abspath(l[1])
        elif 'annotation' in l[0]:
            d_args['annotation'] = os.path.abspath(l[1])
        elif 'parent' in l[0]:
            d_args['parent'] = os.path.abspath(l[1])
        elif l[0].startswith('#'):
            continue
        else:
            #d_args['some_experiment_name'] = [condition, sample_id]
            d_args[l[0]].append([l[2], l[1]])
    return d_args


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group('required arguments')
    required.add_argument("-c", "--config", dest="config",
            help="Configuration file that was passed to Cuffmerge")
    args = p.parse_args()
    if not (args.config):
        p.print_help()
        sys.exit(1)
    cuffdiff(args.config)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()