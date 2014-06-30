#!/usr/bin/env bash
#BSUB -J trim[1-20]
#BSUB -e trim.%J.%I.err
#BSUB -o trim.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P duval_canine

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/canine/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fqin=$DATA/$sample.fastq.gz
fqout=$DATA/$sample.trm.fastq.gz

fastq-mcf $ADAPTORS $fqin -o $fqout -l 15 -L 25 -k 0
