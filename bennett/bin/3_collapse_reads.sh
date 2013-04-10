#! /usr/bin/env bash
#BSUB -J collapse_reads[1-6]
#BSUB -e collapse_reads.%J.%I.err
#BSUB -o collapse_reads.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

<<DOC
Collapse reads into their unique subsets.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
results=$RESULTS/$sample
joined=$results/$sample.joined.fastq.gz
collapsed=$results/$sample.collapsed.fastq.gz
if [[ ! -f $collapsed ]]; then
    python $BIN/collapse_reads.py $joined | gzip -c > $collapsed
fi
