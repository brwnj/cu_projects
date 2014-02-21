#!/usr/bin/env bash
#BSUB -J cdr3_local_alignment
#BSUB -e cdr3_local_alignment.%J.err
#BSUB -o cdr3_local_alignment.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

python $HOME/projects/bennett/bin/target_alignment.py -l 2 -i .50 -q $RESULTS/5_*.txt -t $RESULTS/peptide_seqs.tsv > $RESULTS/alignments.m50.l2.tsv
