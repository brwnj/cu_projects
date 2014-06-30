#!/usr/bin/env bash
#BSUB -J bam2fastq[1-24]
#BSUB -e bam2fastq.%J.%I.err
#BSUB -o bam2fastq.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1

<<DOC
convert alignments to TB from bam to fastq in order to realign to hg19.
DOC

set -o nounset -o pipefail -o errexit -x

samples=(E1T1_Inf E1T1_Uninf E1T24_Inf E1T24_Uninf E1T2_Inf E1T2_Uninf
            E1T8_Inf E1T8_Uninf E2T1_Inf E2T1_Uninf E2T24_Inf E2T24_Uninf
            E2T2_Inf E2T2_Uninf E2T8_Inf E2T8_Uninf E3T1_Inf E3T1_Uninf
            E3T24_Inf E3T24_Uninf E3T2_Inf E3T2_Uninf E3T8_Inf E3T8_Uninf)
sample=${samples[$(($LSB_JOBINDEX - 1))]}

data=$HOME/projects/walter/data/20130204
if [ ! -d $data ]; then
    mkdir -p $data
fi
fastq=$data/$sample.tb.fq
bam=$HOME/projects/walter/results/common/$sample/$sample.tb.novo.bam

java -jar ~/opt/picard-tools-1.79/SamToFastq.jar I=$bam F=$fastq RC=true

gzip $fastq