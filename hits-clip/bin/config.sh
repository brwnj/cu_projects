#!/usr/bin/env bash
set -o nounset

SAMPLES=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
            MP38 MP39ACTG MP39TCGA MP40 MP41 MP42ACTG MP42TCGA MP43ACTG
            MP43TCGA MP44ACTG MP44TCGA MP45ACTG MP45TCGA MP7 MP9 PK11
            PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53
            PK54 helaa helab PK61 PK62 PK63 MCF10A PK12-0 PK12-24 PK12-6)

PI=pillai_kabos_hitsclip
PROJECT=hits-clip
RESULTS=$HOME/projects/hits-clip/results/common
DATA=$HOME/projects/hits-clip/data/common
SRC=$HOME/projects/hits-clip/bin/scripts
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
FASTAS=$HOME/ref/hg19/fa_per_chr
SIZES=$HOME/ref/hg19/hg19.sizes
GENOMEDATA=$RESULTS/hitsclip.gd
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
AGOFOOTPRINT=60
HUB=$RESULTS/hub
GENOME=hg19

SAMPLE_GROUPS=(
BMEC
BT474
BT474estr
BT474herc
HELA
hMSC
HS27A
HS5
HUVEC
MCF7
MCF7estr
MDA231
MCF7_E_72
MCF7_TAM_72
MCF7_TAMRES
MCF10A_group
PK12-0_group
PK12-6_group
PK12-24_group
)
GROUP_REPLICATES=(
"MP42ACTG MP45ACTG MP45TCGA"
"PK21 PK41 PK52"
"PK22 PK53"
"PK23"
"helaa helab"
"MP36 MP43ACTG MP43TCGA MP44ACTG MP44TCGA"
"MP1 MP21 MP35"
"MP2 MP20 MP34"
"MP24 MP38"
"PK11 PK31 PK51"
"PK12 PK32"
"PK24 PK42 PK54"
"PK61"
"PK62"
"PK63"
"MCF10A"
)

if [[ ! -d $HUB ]]; then
    mkdir -p $HUB
fi