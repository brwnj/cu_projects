#! /usr/bin/env bash
#BSUB -J dexseq_main[1-2]
#BSUB -e dexseq_main.%J.%I.err
#BSUB -o dexseq_main.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

samples=(MP PK)
sample=${samples[$(($LSB_JOBINDEX - 1))]}
counts=$RESULT/${sample}*/*.counts.txt.gz

cd $DEXRESULTS

# run_dexseq.py submits jobs directly into the queue
python $BIN/run_dexseq.py -p pillai_kabos_polya -q normal $RUNDEXSEQ $counts
