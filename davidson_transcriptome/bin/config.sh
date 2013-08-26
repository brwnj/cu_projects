#!/usr/bin/env bash
set -o nounset

# 10
SAMPLES=(1_control 1_cytokine 2_control 2_cytokine 3_control 3_cytokine 4_control 4_cytokine 5_control 5_cytokine)
PI=davidson
RESULTS=$HOME/projects/davidson_transcriptome/results/common
DATA=$HOME/projects/davidson_transcriptome/data/common
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
GENES=$HOME/ref/hg19/hg19.gtf
SEQUENCE=$HOME/ref/hg19/hg19.fa
