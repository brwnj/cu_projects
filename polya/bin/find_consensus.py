#!/usr/bin/env python
# encoding: utf-8
"""
Find the consensus of peaks among samples.

Assumes sample name is first in the file name, delimited by either "." or "_"
from the rest of the file name.
"""
import os, sys, tempfile
import os.path as op
import subprocess as sp
from toolshed import reader

def filter_peaks(files, classes):
    """removes peaks that do not match a class in classes."""
    tmps = {}
    classes = set(classes)
    for f in files:
        sample = op.basename(f).split(".")[0].split("_")[0]
        tmp = open(tempfile.mkstemp(suffix=".bed")[1], 'w')
        res = ["chrom", "start", "stop", "name", "score", "strand"]
        for l in reader(f, header=res):
            c = int(l['name'].split(":")[1])
            if c not in classes: continue
            tmp.write("\t".join(l[i] for i in res) + "\n")
        tmp.close()
        tmps[sample] = tmp.name
    return tmps

def add_slop(files, sizes, n):
    """add slop onto bed regions for each peak file.
    files = {sample_name:file_path}
    """
    tmps = {}
    for sample, f in files.iteritems():
        tmp = tempfile.mkstemp(suffix=".bed")[1]
        cmd = "bedtools slop -b %d -i %s -g %s > %s" % (n, f, sizes, tmp)
        sp.call(cmd, shell=True)
        tmps[sample] = tmp
    return tmps

def cleanup(files):
    for f in files.values():
        os.remove(f)

def intersect(files, fraction):
    """files = {sample_name:file_path}"""
    fraction = int(len(files) * fraction)
    names = " ".join(files.keys())
    paths = " ".join(files.values())
    cmd = "|bedtools multiinter -cluster -header -names %s -i %s" % (names, paths)
    for l in reader(cmd, header=True):
        if int(l['num']) < fraction: continue
        print "\t".join([l['chrom'], l['start'], l['end']])

def main(args):
    tmps = filter_peaks(args.files, args.classes)
    tmps_with_slop = add_slop(tmps, args.sizes, args.bases)
    cleanup(tmps)
    intersect(tmps_with_slop, args.fraction)
    cleanup(tmps_with_slop)

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__,
            formatter_class=argparse.RawDescriptionHelpFormatter)

    p.add_argument("sizes", help="chromosome sizes for specific genome")
    p.add_argument("files", nargs="+", help="classified peaks")

    psites = p.add_argument_group("poly(A) sites")
    psites.add_argument("-b", dest="bases", type=int, default=5,
            help="increase region -b base pairs in each direction [%(default)s]")
    psites.add_argument("-c", metavar="CLASS", dest="classes", action="append",
            type=int, default=[1], choices=[1,2,3,4],
            help="class of peaks used to generate consensus %(default)s")

    pinter = p.add_argument_group("intersecting")
    pinter.add_argument("-f", dest="fraction", type=float, default=0.25,
            help="fraction of samples where peak must have been called [%(default)s]")

    args = p.parse_args()
    main(args)
