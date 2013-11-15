#!/usr/bin/env bash
#BSUB -J freebayes[1-2]
#BSUB -e freebayes.%J.%I.err
#BSUB -o freebayes.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P hanson_variants

<<DOC
call variants using freebayes
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/hanson_variants/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

bam=$RESULTS/$sample/$sample.bam
vcf=${bam/.bam/.vcf}
if [[ ! -f $vcf ]]; then
    freebayes -b $bam -v $vcf -f $MM9_FASTA -0
fi
