#!/usr/bin/env python
# encoding: utf-8
"""
Creates a manifest and runs Cuffmerge per condition defined in config.txt. The
config file can contain multiple experiments but each line of a particular
experiment has to start with the same experiment name.

Created by Joe Brown on 2011-12-16.
"""
import argparse
import os
import sys
import gzip
import bz2
from itertools import izip
import urllib
import fnmatch
from subprocess import Popen
from collections import defaultdict


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


def create_manifest(manifest, parent, samples):
    """Writes the manifest of the transcripts.gtf files to include in
    Cuffmerge.
    """
    sh = open(manifest, 'w')
    file_dirs = []
    for root, dirs, files in os.walk(parent):
        for f in fnmatch.filter(files, 'transcripts.gtf'):
            filepath = os.path.join(root, f)
            for s in samples:
                if s in filepath:
                    file_dirs.append(filepath)
    [[sh.write('%s\n' % d)] for d in file_dirs]
    sh.close()
    return True


def run_cuffmerge(config):
    """Creates the manifest and runs Cuffmerge."""
    for k, v in config.iteritems():
        if not k in ('index','annotation','parent'):
            samples = []
            for sample in v:
                samples.append(sample[1])
            man_path = '%s/%s' % (config['parent'], k)
            if not os.path.exists(man_path):
                os.mkdir(man_path)
            manifest = '%s/manifest.txt' % man_path
            cuffmerge = '%s/cuffmerge.sh' % man_path
            if create_manifest(manifest, config['parent'], samples):
                sh = open(cuffmerge, 'w')
                sh.write('#!/bin/sh\n')
                sh.write('#PBS -l walltime=48:00:00,nodes=1:ppn=8,mem=30gb\n')
                sh.write('#PBS -j oe\n')
                sh.write('cd %s\n' % man_path)
                sh.write('cuffmerge --ref-gtf %s --num-threads 8 --ref-sequence %s %s'\
                    % (config['annotation'], config['index'], manifest))
                sh.close()
                Popen(['qsub', cuffmerge])
            else:
                sys.stderr.write("Manifest failed to write while processing: %s\n" % k)
                sys.exit(1)


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


def write_config():
    """Writes a sample configuration file."""
    conf = open('%s/sample_configuration.txt' % os.getcwd(), 'w')
    conf.write("#Comment lines are ignored.\n")
    conf.write("#Keep variable names the same.\n")
    conf.write("#Multiple experiments can be listed.\n")
    conf.write("index\t/mnt/storage/reference/genome.fa\n")
    conf.write("annotation\t/mnt/storage/reference/genes.gtf\n")
    conf.write("parent\t/mnt/storage/experiment/cufflinks/\n")
    conf.write("cancer_control\t012345\tcancer\n")
    conf.write("cancer_control\t012346\tcancer\n")
    conf.write("cancer_control\t987654\tcontrol\n")
    conf.write("cancer_control\t987653\tcontrol\n")


def main():
    p = argparse.ArgumentParser(description=__doc__,
                   formatter_class=argparse.RawDescriptionHelpFormatter)
    required = p.add_argument_group('required arguments')
    required.add_argument("-c", dest="config", 
                            help="File containing experiment specifics")
    p.add_argument("--sample", dest="sample", action="store_true", 
                            help="Write a sample config file")
    args = p.parse_args()
    if args.config:
        config = get_args(args.config)
        run_cuffmerge(config)
    elif args.sample:
        write_config()
    else:
        p.print_help()
        sys.exit(1)


if __name__ == "__main__":
    import doctest
    if doctest.testmod(optionflags=doctest.ELLIPSIS |\
                                   doctest.NORMALIZE_WHITESPACE).failed == 0:
        main()