#!/usr/bin/env bash
#BSUB -J maketagdir[1-20]
#BSUB -e maketagdir.%J.%I.err
#BSUB -o maketagdir.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

<<DOC
write tags into individual result folders
DOC

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
bam=$RESULTS/$sample/$sample.bam
tagdir=$RESULTS/$sample

makeTagDirectory $tagdir $bam -genome hg19 -normGC default
