#!/usr/bin/env bash
set -o nounset

SAMPLES=(PK11 PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51
PK52 PK53 PK54 helaa helab PK61 PK62 PK63 MCF10A PK12-0 PK12-24
PK12-6
MCF7-0
MCF7-6
MCF7-24
MCF7-0_CAC
MCF7-24_GAT
MCF7-6_ATT
MCF7-0_TGT
MCF7-24_CTG
MCF7-6_GGC
)

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
HUBNAME=HITSCLIP_201405

TRIMSCRIPT=$SRC/trim_adapter.py
TRIMADAPTER=GTGTCA

FILTERSCRIPT=$SRC/filter_bam.py
FILTERDOWNSTREAMBASES=25

PEAKSNONUNIQUEQ=0.00001
PEAKSUNIQUEQ=0.01

SEEDFASTA8=/vol1/home/brownj/ref/mirbase/20/mature.hsa.8mer_seed.fa.gz
SEEDFASTA7=/vol1/home/brownj/ref/mirbase/20/mature.hsa.7mer_seed.fa.gz
GENEANNOTATION=/vol1/home/brownj/ref/hg19/refseq.wholegene.bed.gz
MIRBASEREGIONS=/vol1/home/brownj/ref/mirbase/20/mature.hsa.bed.gz

declare -A REPLICATES
REPLICATES=([BT474]="PK21 PK41 PK52"
		[BT474estr]="PK22 PK53"
		[HELA]="helaa helab"
		[MCF7]="PK11 PK31 PK51"
		[MCF7estr]="PK12 PK32"
		[MDA231]="PK24 PK42 PK54"
)

# lone samples or i just don't know which group they belong to
SINGLETONS=(PK23
PK33
PK61
PK62
PK63
MCF10A
PK12-0
PK12-24
PK12-6
MCF7-0
MCF7-6
MCF7-24
MCF7-0_CAC
MCF7-24_GAT
MCF7-6_ATT
MCF7-0_TGT
MCF7-24_CTG
MCF7-6_GGC
)

# singletons as well as group names
MAPSEEDS=(PK23
PK33
PK61
PK62
PK63
MCF10A
PK12-0
PK12-24
PK12-6
MCF7-0
MCF7-6
MCF7-24
MCF7-0_CAC
MCF7-24_GAT
MCF7-6_ATT
MCF7-0_TGT
MCF7-24_CTG
MCF7-6_GGC
BT474
BT474estr
HELA
MCF7
MCF7estr
MDA231
)

declare -A COLORS
COLORS=([PK11]="0,109,44"
        [PK31]="35,139,69"
        [PK51]="65,174,118"
        [MCF7]="0,68,27"
        [PK21]="8,81,156"
        [PK41]="33,113,181"
        [PK52]="66,146,198"
        [BT474]="8,48,107"
        [PK24]="166,54,3"
        [PK42]="217,72,1"
        [PK54]="241,105,19"
        [MDA231]="127,39,4"
		[PK12]="63,0,125"
		[PK22]="8,29,88"
		[PK23]="73,0,106"
		[PK32]="221,52,151"
		[PK33]="239,59,44"
		[PK53]="188,189,220"
		[helaa]="103,169,207"
		[helab]="54,144,192"
		[PK61]="254,178,76"
		[PK62]="199,233,180"
		[PK63]="247,129,191"
		[MCF10A]="122,1,119"
		[PK12-0]="204,76,2"
		[PK12-24]="153,52,4"
		[PK12-6]="102,37,6"
		[MCF7-0]="231,41,138"
		[MCF7-6]="173,221,142"
		[MCF7-24]="128,125,186"
		[MCF7-0_CAC]="231,41,138"
		[MCF7-24_GAT]="128,125,186"
		[MCF7-6_ATT]="173,221,142"
		[MCF7-0_TGT]="231,41,138"
		[MCF7-24_CTG]="128,125,186"
		[MCF7-6_GGC]="173,221,142"
)
