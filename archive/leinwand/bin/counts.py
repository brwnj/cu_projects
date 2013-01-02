#! /usr/bin/env python
# encoding: utf-8
from bsub import bsub
from pybedtools import BedTool
from toolshed import reader
import os
import os.path as op
import pandas as pd
import sys
import fnmatch

def getfilelist(path, pattern):
    files = []
    for root, dirnames, filenames in os.walk(path):
          for filename in fnmatch.filter(filenames, pattern):
              files.append(op.join(root, filename))
    return files

def counts(samples, result_path, peak_ext, bam_ext):
    # get the consensus peaks
    f = open("%s/peak_coordinates.bed" % result_path, 'w')
    x = BedTool()
    consensus = x.multi_intersect(i=getfilelist(result_path, "*%s" % peak_ext))
    for c in consensus:
        # fixing formatting from bedtool object
        replicate_counts = c.name
        if replicate_counts < 2: continue
        
        fields = [c.chrom, c.start, c.stop, "%s:%d-%d\n" % \
                    (c.chrom, c.start, c.stop)]
        f.write("\t".join(map(str, fields)))
    f.close()
    # get counts for each sample
    jobs = []
    countfiles = []
    for sample in samples:
        bams = getfilelist(result_path, sample + "*%s" % bam_ext)
        assert(len(bams) == 1)
        outdir = result_path.rstrip("/") + "/" + sample
        countsresult = outdir + "/" + sample + ".counts"
        countfiles.append(countsresult)
        if op.exists(countsresult): continue
        cmd = "bedtools coverage -abam %s -b %s > %s" % \
                    (bams[0], f.name, countsresult)
        jobid = bsub(sample + "_counts", 
                        R="select[mem>16] rusage[mem=16] span[hosts=1]",
                        verbose=True)(cmd)
        jobs.append(jobid)
    bsub.poll(jobs)
    # counts to matrix
    allcounts = {}
    for cf in countfiles:
        cfname = op.basename(cf).split(".counts")[0]
        casecounts = {}
        for toks in reader(cf, header="chrom start stop name a_overlaps_in_b \
                    b_with_nonzero length_b frac_b_nonzero".split()):
            casecounts[toks['name']] = int(toks['a_overlaps_in_b'])
        allcounts[cfname] = casecounts
    countsdf = pd.DataFrame(allcounts)
    countsdf.to_csv(sys.stdout, sep="\t", header=True)

def main(args):
    counts(args.samples, args.result_path, args.peak_ext, args.alignment_ext)

if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--samples', required=True, nargs="+", help='sample names')
    p.add_argument('--result_path', required=True, help='parent directory of sample results')
    p.add_argument('--peak_ext', default="_peaks.bed", help="file extension to search for [ _peaks.bed ]")
    p.add_argument('--alignment_ext', default=".bam", help="file extension to search for [ .bam ]")
    args = p.parse_args()
    main(args)