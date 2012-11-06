#! /usr/bin/env bash
#BSUB -J fastqc[1-7]
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

<<DOC
qc the reads
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(workaround
1_TCAG_L007_R1_001.fastq.gz
2_AGCA_L007_R1_001.fastq.gz
3_CTGT_L007_R1_001.fastq.gz
4_GATC_L007_R1_001.fastq.gz
5_TGAG_L007_R1_001.fastq.gz
6_ACTC_L007_R1_001.fastq.gz
7_CAGT_L007_R1_001.fastq.gz
)

SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}

PROJECT=$HOME/projects/huang
OUT=$PROJECT/results/common/$(basename $SAMPLE .fastq.gz)

mkdir -p $OUT

$HOME/opt/FastQC/fastqc --outdir=$OUT \
    $PROJECT/data/20120615/$SAMPLE