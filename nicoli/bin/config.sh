#!/usr/bin/env bash
set -o nounset

SAMPLES=(L1 L2)
RESULTS=$HOME/projects/nicoli/results/common
DATA=$HOME/projects/nicoli/data/20131001
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
FASTAS=$HOME/ref/hg19/fa_per_chr
SIZES=$HOME/ref/hg19/hg19.sizes
GENOMEDATA=$RESULTS/nicoli.gd
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
HUB=$RESULTS/hub
