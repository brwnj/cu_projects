#!/usr/bin/env bash
#BSUB -J sort_umi[1-18]
#BSUB -e sort_umi.%J.%I.err
#BSUB -o sort_umi.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1=$READS/${sample}_R1.umi.fastq.gz
r2=$READS/${sample}_R2.umi.fastq.gz

r1_out=$READS/${sample}_R1.umi.sorted.fastq.gz
r2_out=$READS/${sample}_R2.umi.sorted.fastq.gz

umi_length=8
bin=$HOME/projects/bennett/bin

if [[ ! -f $r1_out ]]; then
    python $bin/process_umi.py sort $r1 $umi_length | gzip -c > $r1_out
fi
if [[ ! -f $r2_out ]]; then
    python $bin/process_umi.py sort $r2 $umi_length | gzip -c > $r2_out
fi
