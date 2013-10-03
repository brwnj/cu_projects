#!/usr/bin/env bash
#BSUB -J reproc.bw[1-6]
#BSUB -e reproc.bw.%J.%I.err
#BSUB -o reproc.bw.%J.%I.out
#BSUB -q normal

<<DOC
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/polya/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULT/$sample

umibam=$RESULT/$sample/$sample.UMIs_not_removed.bam
bam=$results/$sample.bam

bam2bw.py -5 -b -v $umibam $CHROM_SIZES $PROJECTID
bam2bw.py -5 -b -v $bam $CHROM_SIZES $PROJECTID
