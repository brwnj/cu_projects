#! /usr/bin/env bash
#BSUB -J fastqc[1-20]
#BSUB -e fastqc.%J.%I.err
#BSUB -o fastqc.%J.%I.out
#BSUB -q normal

<<DOC
qc the reads
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(workaround
10_AAGCTA_L006_R1_001.fastq.gz
11_GTAGCC_L006_R1_001.fastq.gz
12_TACAAG_L006_R1_001.fastq.gz
13_TTGACT_L006_R1_001.fastq.gz
14_GGAACT_L006_R1_001.fastq.gz
15_TGACAT_L006_R1_001.fastq.gz
16_GGACGG_L006_R1_001.fastq.gz
17_GCGGAC_L006_R1_001.fastq.gz
18_TTTCAC_L006_R1_001.fastq.gz
19_GGCCAC_L006_R1_001.fastq.gz
1_CGTGAT_L006_R1_001.fastq.gz
20_CGAAAC_L006_R1_001.fastq.gz
21_TGGTCA_L006_R1_001.fastq.gz
22_TCAAGT_L006_R1_001.fastq.gz
23_ACATCG_L006_R1_001.fastq.gz
24_ATTGGC_L006_R1_001.fastq.gz
3_GCCTAA_L006_R1_001.fastq.gz
5_CACTGT_L006_R1_001.fastq.gz
7_GATCTG_L006_R1_001.fastq.gz
9_CTGATC_L006_R1_001.fastq.gz
)

SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}

PROJECT=$HOME/projects/thaikoottathil
OUT=$PROJECT/results/common/$(basename $SAMPLE .fastq.gz)

mkdir -p $OUT

$HOME/opt/FastQC/fastqc --outdir=$OUT \
    $PROJECT/data/$SAMPLE