#!/usr/bin/env bash
set -o nounset

SAMPLES=(L1_001 L1_002 L1_003 L1_004 L1_005 L1_006 L1_007 L2_001 L2_002 L2_003 L2_004 L2_005 L2_006 L2_007)
RESULTS=$HOME/projects/nicoli/results/common
DATA=$HOME/projects/nicoli/data/20131001
NOVOIDX=$HOME/projects/hits-clip/data/common/novoalign/hg18
GENOMEDATA=$RESULTS/nicoli.rmd.gd
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
SIZES=$HOME/ref/hg18/hg18.sizes
FASTAS=$HOME/ref/hg18/fa_per_chr
HUB=$RESULTS/hub
