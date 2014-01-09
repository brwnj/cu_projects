#!/usr/bin/env bash
#BSUB -J findpeaks[1-3]
#BSUB -e findpeaks.%J.%I.err
#BSUB -o findpeaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P cohrs_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/cohrs_chipseq/bin/config.sh

groups=(H3K9ac H3PAN R_IgG)
group=${groups[$(($LSB_JOBINDEX - 1))]}

findPeaks $RESULTS/$group -style factor -o auto -i $RESULTS/INPUT
