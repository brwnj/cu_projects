#!/usr/bin/env python
# encoding: utf-8
"""
Leinwand
"""
from bsub import bsub
import sys

sys.path.append('/vol1/home/brownj/projects/utils')
import ngseq

def main(args):
    samples = ["MDX_22_AGTTCC_L003_R1_001",
               "MDX_23_ATGTCA_L003_R1_001",
               "MDX_24_CCGTCC_L003_R1_001",
               "WT_21_AGTCAA_L003_R1_001",
               "WT_25_GTAGAG_L003_R1_001",
               "WT_42_GTCCGC_L003_R1_001"]
    datadir = "/vol1/home/brownj/projects/leinwand/data/20121101"
    adapters = "%s/adapters.fa" % datadir
    resultsdir = "/vol1/home/brownj/projects/leinwand/results/common"
    fastqc_script = "/vol1/home/brownj/opt/fastqc/fastqc"
    picard = "/vol1/home/brownj/opt/picard-tools-1.79"
    reference_fasta = "/vol1/home/brownj/ref/mm9/mm9.fa"
    gmapdb = "/vol1/home/brownj/ref/gmapdb"
    gsnapcmd = "gsnap -D {} -d mm9 --gunzip \
                --batch=5 --nofails --nthreads=4 --format=sam -v snp128_strict_wholeChrs {} \
                | samtools view -ShuF 4 - \
                | samtools sort -o - {}.temp -m 9500000000 > {}"
    chrom_sizes = "/vol1/home/brownj/ref/mm9/mm9.sizes"
    
    if args.clobber:
        ngseq.clobber_previous(resultsdir)
    
    ngseq.fastqc(fastqc_script, samples, datadir)
    bsub.poll(ngseq.trimadapter(datadir, adapters))
    bsub.poll(ngseq.gsnap(samples, datadir, resultsdir, gmapdb, gsnapcmd))
    # ngseq.alignment_stats(resultsdir, picard, reference_fasta)
    ngseq.cleanup(resultsdir)

    # create genomedata archive in results/common
    bam_pattern = "/vol1/home/brownj/projects/leinwand/results/common/*/*.bam"
    output_dir = "/vol1/home/brownj/projects/leinwand/results/common"
    # ngseq.bam2gd(chrom_sizes, reference_fasta, bam_pattern, output_dir)
    
    # call peaks using genomedata archive

    # combine positive and negative strand
    # 
    # ngseq.counts(samples, result_path, peak_ext="_peaks.bed", bam_ext="bam")

if __name__ == '__main__':
    import argparse
    p = argparse.ArgumentParser(description=__doc__, 
                        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--clobber", action="store_true", 
                        help="clear all previous results")
    args = p.parse_args()
    main(args)