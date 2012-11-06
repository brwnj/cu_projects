#! /usr/bin/env bash
#BSUB -J peaks[1-6]
#BSUB -e peaks.%J.%I.err
#BSUB -o peaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1

<<DOC
use macs to call peaks.
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(workaround
2_AGCA_L007_R1_001.bam
3_CTGT_L007_R1_001.bam
4_GATC_L007_R1_001.bam
5_TGAG_L007_R1_001.bam
6_ACTC_L007_R1_001.bam
7_CAGT_L007_R1_001.bam
)

SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}
NAME=$(basename $SAMPLE .bam)
PROJECT=$HOME/projects/huang

CONTROL=$PROJECT/results/common/1_TCAG_L007_R1_001/1_TCAG_L007_R1_001.bam

#if wiggle dir exists, this will fail
macs14 --treatment $PROJECT/results/common/$(basename $SAMPLE .bam)/$SAMPLE \
    --control $CONTROL \
    --name $NAME \
    --format BAM \
    --gsize 12000000 \
    --pvalue 1e-3
    --wig \
    --single-profile \
    --call-subpeaks