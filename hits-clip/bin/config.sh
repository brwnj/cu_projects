#!/usr/bin/env bash
set -o nounset

SAMPLEIDS=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
            MP38 MP39ACTG MP39TCGA MP40 MP41 MP42ACTG MP42TCGA MP43ACTG
            MP43TCGA MP44ACTG MP44TCGA MP45ACTG MP45TCGA MP7 MP9 PK11
            PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53
            PK54 helaa helab PK61 PK62 PK63)

PI=pillai_kabos_hitsclip
PROJECT=hits-clip
RESULTS=$HOME/projects/hits-clip/results/common
DATA=$HOME/projects/hits-clip/data/common
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
FASTAS=$HOME/ref/hg19/fa_per_chr
SIZES=$HOME/ref/hg19/hg19.sizes
GENOMEDATA=$RESULTS/hitsclip.gd
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)