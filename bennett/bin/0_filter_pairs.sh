#! /usr/bin/env bash
#BSUB -J filter_pairs[1-6]
#BSUB -e filter_pairs.%J.%I.err
#BSUB -o filter_pairs.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

<<DOC
filter out paired-end reads that have no mate.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
r1=$READS/${sample}_R1.fastq.gz
r2=$READS/${sample}_R2.fastq.gz
r1filtered=$READS/${sample}_R1.filtered.fastq
r2filtered=$READS/${sample}_R2.filtered.fastq

if [[ ! -f $r1filtered ]]; then
    # output written to same dir as r1 and r2; includes "filtered" in name
    python $BIN/filter_pairs.py $r1 $r2
    gzip $r1filtered $r2filtered
fi