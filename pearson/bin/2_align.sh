#!/usr/bin/env bash
#BSUB -J novo[1-2]
#BSUB -e novo.%J.%I.err
#BSUB -o novo.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>21] rusage[mem=20] span[hosts=1]"
#BSUB -n 8
#BSUB -P pearson

<<DOC
align paired-end samples using novoalign
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/pearson/bin/config.sh

pair=${PAIRS[$(($LSB_JOBINDEX - 1))]}
name=${NAMES[$(($LSB_JOBINDEX - 1))]}

outdir=$RESULTS/common/$name
bam=$outdir/$name.remac.bam
stats=$outdir/$name.remac.alignment_stats.txt
idx=$HOME/ref/tetrahymena/JCVI_TTA1_2_2_18.novoidx

if [ ! -f $bam ]; then
    novoalign -d $idx -f $pair -o SAM -r None -i 250 100 -c 8 -k \
        2> $stats \
        | samtools view -ShuF4 - \
        | samtools sort -o - $name.temp -m 8G \
        > $bam
fi