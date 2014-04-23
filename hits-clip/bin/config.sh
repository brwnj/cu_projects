#!/usr/bin/env bash
set -o nounset

SAMPLES=(PK11 PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51
PK52 PK53 PK54 helaa helab PK61 PK62 PK63 MCF10A PK12-0 PK12-24
PK12-6 MCF7-0 MCF7-6 MCF7-24)

PI=hitsclip
PROJECT=hits-clip
RESULTS=$HOME/projects/hits-clip/results/common
DATA=$HOME/projects/hits-clip/data/common
SRC=$HOME/projects/hits-clip/bin/scripts
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
FASTAS=$HOME/ref/hg19/fa_per_chr
FASTA=$HOME/ref/hg19/hg19.fa
SIZES=$HOME/ref/hg19/hg19.sizes
GENOMEDATA=$RESULTS/hitsclip.gd
AGOFOOTPRINT=60
HUB=$RESULTS/hub
GENOME=hg19

TRIMSCRIPT=/vol1/home/brownj/projects/hits-clip/bin/scripts/trim_adapter.py
TRIMADAPTER=GTGTCA

FILTERSCRIPT=/vol1/home/brownj/projects/hits-clip/bin/scripts/filter_bam.py
FILTERDOWNSTREAMBASES=25

PEAKSNONUNIQUEQ=0.00001
PEAKSUNIQUEQ=0.01

declare -A REPLICATES
REPLICATES=([BT474]="PK21 PK41 PK52"
		[BT474estr]="PK22 PK53"
		[BT474herc]="PK23"
		[HELA]="helaa helab"
		[MCF7]="PK11 PK31 PK51"
		[MCF7estr]="PK12 PK32"
		[MDA231]="PK24 PK42 PK54"
		[MCF7_E_72]="PK61"
		[MCF7_TAM_72]="PK62"
		[MCF7_TAMRES]="PK63"
		[MCF10A_group]="MCF10A"
		[PK12-0_group]="PK12-0"
		[PK12-6_group]="PK12-6"
		[PK12-24_group]="PK12-24")
