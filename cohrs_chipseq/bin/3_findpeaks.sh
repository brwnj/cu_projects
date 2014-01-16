#!/usr/bin/env bash
#BSUB -J findpeaks
#BSUB -e findpeaks.%J.err
#BSUB -o findpeaks.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P cohrs_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/cohrs_chipseq/bin/config.sh

# groups=(H3K9ac H3PAN R_IgG)
# group=${groups[$(($LSB_JOBINDEX - 1))]}

findPeaks $RESULTS/H3K9ac -style factor -o auto -i $RESULTS/INPUT
findPeaks $RESULTS/H3PAN -style histone -o auto -i $RESULTS/INPUT
findPeaks $RESULTS/R_IgG -style histone -o auto -i $RESULTS/INPUT
