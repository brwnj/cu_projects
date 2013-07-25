#!/usr/bin/env bash
#BSUB -J join_reads[1-18]
#BSUB -e join_reads.%J.%I.err
#BSUB -o join_reads.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1=$READS/${sample}_R1.umi.sorted.umifiltered.rmadptr.fastq.gz
r2=$READS/${sample}_R2.umi.sorted.umifiltered.rmadptr.fastq.gz
joined=$READS/${sample}.joined.fastq.gz

bin=$HOME/projects/bennett/bin

if [[ ! -f $joined ]]; then
    python $bin/join_reads.py $r1 $r2 | gzip -c > $joined
fi
