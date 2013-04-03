#! /usr/bin/env bash
#BSUB -J trim_primer[1-6]
#BSUB -e trim_primer.%J.%I.err
#BSUB -o trim_primer.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

<<DOC
Trim the primer sequence from the 5' and 3' ends of read 1 and read 2, adding
the primer name to the read name.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
r1=$READS/${sample}_R1.filtered.fastq.gz
r2=$READS/${sample}_R2.filtered.fastq.gz
r1trimmed=$READS/${sample}_R1.filtered.trimmed.fastq
r2trimmed=$READS/${sample}_R2.filtered.trimmed.fastq

if [[ ! -f $r1trimmed.gz ]]; then
    python $BIN/trim_adapters.py $r1 $r2 $R1PRIMERS $R2PRIMERS $r1trimmed $r2trimmed
    gzip $r1trimmed $r2trimmed
fi