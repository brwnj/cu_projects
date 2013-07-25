#!/usr/bin/env bash
#BSUB -J trim_adapters[1-18]
#BSUB -e trim_adapters.%J.%I.err
#BSUB -o trim_adapters.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1=$READS/${sample}_R1.umi.sorted.umifiltered.fastq.gz
r2=$READS/${sample}_R2.umi.sorted.umifiltered.fastq.gz

r1_out=$READS/${sample}_R1.umi.sorted.umifiltered.rmadptr.fastq.gz
r2_out=$READS/${sample}_R2.umi.sorted.umifiltered.rmadptr.fastq.gz

if [[ ! -f $r1_out ]] || [[ ! -f $r2_out ]]; then
    python $BIN/trim_adapters.py $r1 $r2 $R1PRIMERS $R2PRIMERS
fi
