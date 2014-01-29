#!/usr/bin/env bash
#BSUB -J cdr3_local_alignment[1-5]
#BSUB -e cdr3_local_alignment.%J.%I.err
#BSUB -o cdr3_local_alignment.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

mismatches=(".85" ".75" ".65" ".55" ".50")
mismatch=${mismatches[$(($LSB_JOBINDEX - 1))]}

python $HOME/projects/bennett/bin/cdr3_target_alignment.py -i $mismatch -q $RESULTS/5_*.txt -t $RESULTS/peptide_seqs.tsv > $RESULTS/alignments.$mismatch.tsv
