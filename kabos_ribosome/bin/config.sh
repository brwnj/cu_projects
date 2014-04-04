#!/usr/bin/env bash
PROJECTID=kabos_ribosome
RESULTS=/vol1/home/brownj/projects/kabos_ribosome/results/common
DATA=/vol1/home/brownj/projects/kabos_ribosome/data/common
SAMPLES=(MCF7-0
            MCF7-0_fraction
            MCF7-0_totalRNA
            MCF7-6
            MCF7-6_fraction
            MCF7-6_totalRNA
            PK12-0
            PK12-6)
SIZES=/vol1/home/brownj/ref/hg19/hg19.sizes
TRIMADAPTER=AGATCGGAAGAG
TRIMMINLENGTH=25
FILTERINDEX=/vol1/home/brownj/ref/hg19/rRNA
FILTERSEEDLEN=23
TOPHATREF=/vol1/home/brownj/ref/hg19/hg19
TOPHATSEGMENTMISMATCHES=1
TOPHATTRANSCRIPTOMEINDEX=/vol1/home/brownj/ref/hg19/transcriptome/genes

GENOME=hg19
HUB=$RESULTS/hub
GENOMES=$HUB/genomes.txt
TRACKDB=$HUB/$GENOME/trackDb.txt
HUBNAME=Ribosome_Profiling
HUBLABEL="Ribosome Profiling"
HUBEMAIL=brwnjm@gmail.com

declare -A COLORS
declare -A SUBGROUP1
COLORS=([MCF7-0]="102,194,164"
        [MCF7-0_fraction]="35,139,69"
        [MCF7-0_totalRNA]="0,68,27"
        [MCF7-6]="223,101,176"
        [MCF7-6_fraction]="206,18,86"
        [MCF7-6_totalRNA]="103,0,31"
        [PK12-0]="107,174,214"
        [PK12-6]="253,141,60")
SUBGROUP1=([MCF7-0]=MCF7
            [MCF7-0_fraction]=MCF7
            [MCF7-0_totalRNA]=MCF7
            [MCF7-6]=MCF7
            [MCF7-6_fraction]=MCF7
            [MCF7-6_totalRNA]=MCF7
            [PK12-0]=PK12
            [PK12-6]=PK12)

POSTPROCESSING=$RESULTS/postprocessing
HG19GTF=/vol1/home/brownj/ref/hg19/hg19.gtf
