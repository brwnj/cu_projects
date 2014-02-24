#!/usr/bin/env bash
#BSUB -J "variants[1-24]%8"
#BSUB -e variants.%J.%I.err
#BSUB -o variants.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P leung

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/leung/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

vcf=$RESULTS/$sample/$sample.vcf
snpeff=$RESULTS/$sample/$sample.snpeff.txt.gz

if [[ ! -f $snpeff ]]; then
    java -jar $SNPEFF eff \
        -chr chr \
        -minC 10 \
        -no-downstream \
        -no-intergenic \
        -no-intron \
        -no-upstream \
        -noStats \
        -v \
        -c $SNPEFFCONFIG \
        -o txt \
        hg19 $vcf \
        | gzip -c > $snpeff
fi
