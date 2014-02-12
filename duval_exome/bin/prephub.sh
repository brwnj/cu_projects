#!/usr/bin/env bash
#BSUB -J prephub[1-16]
#BSUB -e prephub.%J.%I.err
#BSUB -o prephub.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P duval_exome


set -o nounset -o pipefail -o errexit -x
source $HOME/projects/duval_exome/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULTS/$sample
bam=$results/$sample.bam
bamindex=$results/$sample.bam.bai
vcf=$results/$sample.vcf
vcfindex=$results/$sample.vcf.idx
bwpos=$results/${sample}_pos.bw

# index the bams
if [[ ! $bamindex ]]; then
    samtools index $bam
fi

# coverage
if [[ ! -f $bwpos ]]; then
    bam2bw.py $bam $SIZES duval_exome
fi

# index the vcfs
if [[ ! -f $vcfindex ]]; then
    ~/opt/IGVTools/igvtools index $vcf
fi
