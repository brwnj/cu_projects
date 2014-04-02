#!/usr/bin/env bash
PROJECTID=kabos_ribosome
RESULTS=/vol1/home/brownj/projects/kabos_ribosome/results/common
DATA=/vol1/home/brownj/projects/kabos_ribosome/data/common
SAMPLES=(
MCF7-0
MCF7-0_fraction
MCF7-0_totalRNA
MCF7-6
MCF7-6_fraction
MCF7-6_totalRNA
PK12-0
PK12-6
)
TRIMADAPTER=AGATCGGAAGAG
TRIMMINLENGTH=25
FILTERINDEX=/vol1/home/brownj/ref/hg19/rRNA
FILTERSEEDLEN=23
TOPHATREF=/vol1/home/brownj/ref/hg19/hg19
TOPHATSEGMENTMISMATCHES=1
TOPHATTRANSCRIPTOMEINDEX=/vol1/home/brownj/ref/hg19/transcriptome/genes

GENOME=hg19
HUB=$RESULTS/hub
declare -A COLORS
declare -A UCSCGROUP
COLORS=([MCF7-0]=""
        [MCF7-0_fraction]=""
        [MCF7-0_totalRNA]=""
        [MCF7-6]=""
        [MCF7-6_fraction]=""
        [MCF7-6_totalRNA]=""
        [PK12-0]=""
        [PK12-6]="")
SUBGROUP1=([PK11]=MCF7
            [PK31]=MCF7
            [PK51]=MCF7
            [PK21]=BT474
            [PK41]=BT474
            [PK52]=BT474
            [PK24]=MDA231
            [PK42]=MDA231
            [PK54]=MDA231
            [MDA231]=MDA231
            [MCF7]=MCF7
            [BT474]=BT474)
