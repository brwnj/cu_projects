#!/usr/bin/env bash
#BSUB -J prepend_umi[1-18]
#BSUB -e prepend_umi.%J.%I.err
#BSUB -o prepend_umi.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

i1=$READS/${sample}_I1.fastq.gz
r1=$READS/${sample}_R1.fastq.gz
r2=$READS/${sample}_R2.fastq.gz

r1_out=$READS/${sample}_R1.umi.fastq.gz
r2_out=$READS/${sample}_R2.umi.fastq.gz

bin=$HOME/projects/bennett/bin

if [[ ! -f $r1_out ]]; then
    python $bin/prepend_umi.py -b 6 -e 14 $i1 $r1 | gzip -c > $r1_out
fi
if [[ ! -f $r2_out ]]; then
    python $bin/prepend_umi.py -b 6 -e 14 $i1 $r2 | gzip -c > $r2_out
fi
