#!/usr/bin/env bash
#BSUB -J trim[1-20]
#BSUB -e trim.%J.%I.err
#BSUB -o trim.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

<<DOC

DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
untrimmed=$DATA/$sample.fastq.gz
trimmed=$DATA/$sample.trimmed.fastq.gz

if [[ ! -f $trimmed ]]; then
    seqtk trimfq $untrimmed | gzip -c > $trimmed
fi
