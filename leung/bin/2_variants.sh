#!/usr/bin/env bash
#BSUB -J variants[1-2]
#BSUB -e variants.%J.%I.err
#BSUB -o variants.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P leung

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/leung/bin/config.sh

SAMPLES=(143_1_20h 143_1HSV1_20h)
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

java='java -Xmx16g -jar'
base_name=$RESULTS/$sample/$sample
bam=$base_name.bam
nodups=$base_name.nodups.bam
duplicatemetrics=$base_name.dup_metrics.txt
vcf=$base_name.vcf

if [[ -f $vcf ]]; then
    echo "processing complete for $sample"
    exit 0
fi

if [ -f $bam ] && [ ! -f $nodups ]; then
    samtools rmdup $bam $nodups
    samtools index $nodups
fi

if [ -f $nodups ] && [ ! -f $vcf ]; then
    freebayes -b $nodups -v $vcf -f $REFERENCE -0
fi
