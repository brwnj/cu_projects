#!/usr/bin/env bash
#BSUB -J concatenate_reads[1-24]
#BSUB -e concatenate_reads.%J.%I.err
#BSUB -o concatenate_reads.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P leung

<<DOC
combine reads across lanes
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/leung/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1=$DATA/common/${sample}_R1.fastq.gz
r2=$DATA/common/${sample}_R2.fastq.gz

zcat $RAWDATA/${sample}*R1*fastq.gz | gzip -c > $r1
zcat $RAWDATA/${sample}*R2*fastq.gz | gzip -c > $r2
