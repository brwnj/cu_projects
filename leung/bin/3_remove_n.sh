#!/usr/bin/env bash
#BSUB -J "qc[1-24]%8"
#BSUB -e qc.%J.%I.err
#BSUB -o qc.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P leung

<<DOC
align rnaseq reads using gsnap
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/leung/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
results=$RESULTS/$sample
orig=$results/$sample.bam
new=$results/$sample.non.bam

samtools view -h $orig | awk '$6!~/N/' | samtools view -Sb - > $new
