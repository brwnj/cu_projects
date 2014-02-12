#!/usr/bin/env bash
#BSUB -J "variants[1-24]%8"
#BSUB -e variants.%J.%I.err
#BSUB -o variants.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P leung

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/leung/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

base_name=$RESULTS/$sample/$sample
bam=$base_name.bam
nonbam=$base_name.non.bam
rmdupbam=$base_name.rmdup.bam
vcf=$base_name.vcf

if [ -f $nonbam ] && [ ! -f $rmdupbam ]; then
    samtools rmdup $nonbam $rmdupbam
    samtools index $rmdupbam
fi

if [ -f $rmdupbam ] && [ ! -f $vcf ]; then
    freebayes -b $rmdupbam -v $vcf -f $REFERENCE -0
fi
