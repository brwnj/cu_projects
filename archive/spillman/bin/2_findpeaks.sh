#!/usr/bin/env bash
#BSUB -J findpeaks[1-2]
#BSUB -e findpeaks.%J.%I.err
#BSUB -o findpeaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P spillman

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/spillman/bin/config.sh

samples=(E21PR_plus E21PR_minus)
controls=(E21PR_plus_input E21PR_minus_input)
sample=${samples[$(($LSB_JOBINDEX - 1))]}
control=${controls[$(($LSB_JOBINDEX - 1))]}

findPeaks $sample -style factor -o auto -i $control
