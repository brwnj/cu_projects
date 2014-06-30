#!/usr/bin/env bash
#BSUB -J transcriptome[1-10]
#BSUB -e transcriptome.%J.%I.err
#BSUB -o transcriptome.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P davidson_transcriptome

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/davidson_transcriptome/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
results=$RESULTS/$sample
bam=$results/$sample.bam

cufflinks -g $GENES -b $SEQUENCE -u -L $sample -o $results -p 8 $bam
