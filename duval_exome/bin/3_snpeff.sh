#!/usr/bin/env bash
#BSUB -J snpeff[1-16]
#BSUB -e snpeff.%J.%I.err
#BSUB -o snpeff.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P duval_exome

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/duval_exome/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

vcf=$RESULTS/$sample/$sample.vcf
snpeff=$RESULTS/$sample/$sample.snpeff.txt.gz
# annotated=$RESULTS/$sample/$sample.annotated_snps.txt.gz

if [[ ! -f $snpeff ]]; then
    java -jar $SNPEFF eff \
        -chr chr \
        -fi $TARGETS \
        -minC 10 \
        -no-downstream \
        -no-intergenic \
        -no-intron \
        -no-upstream \
        -noStats \
        -v \
        -c $SNPEFFCONFIG \
        -o txt \
        canFam3 $vcf \
        | gzip -c > $snpeff
    # python $BIN/annotate.py $snpeff $ANNOTATION | gzip -c > $annotated
fi
