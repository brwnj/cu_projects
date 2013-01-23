#!/usr/bin/env bash
#BSUB -J macs
#BSUB -e macs.%J.err
#BSUB -o macs.%J.out
#BSUB -q normal
#BSUB -R "select[mem>6] rusage[mem=6] span[hosts=1]"
#BSUB -n 1

<<DOC
call peaks using MACS
DOC

results=$HOME/projects/williams/results/common/all
bam=$results/all.bam

macs14 -t $bam -n all -f BAM -g mm -S