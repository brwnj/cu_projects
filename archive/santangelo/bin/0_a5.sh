#!/usr/bin/env bash
#BSUB -J a5[1-7]
#BSUB -o a5.%J.%I.out
#BSUB -e a5.%J.%I.err
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P santangelo

<<DOC
assemble bacterial genome using A5.
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/santangelo/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
r1=$READS/${sample}_R1.fastq.gz
r2=$READS/${sample}_R2.fastq.gz
results=$RESULTS/$sample

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

a5_pipeline.pl $r1 $r2 $sample
