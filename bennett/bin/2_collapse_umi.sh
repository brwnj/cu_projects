#!/usr/bin/env bash
#BSUB -J collapse_umi[1-10]
#BSUB -e collapse_umi.%J.%I.err
#BSUB -o collapse_umi.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1=$READS/${sample}_R1.umi.sorted.fastq.gz
r2=$READS/${sample}_R2.umi.sorted.fastq.gz

umi=NNNNNNNN
bin=$HOME/projects/bennett/bin

python $bin/process_umi.py collapse $r1 $r2 $umi
