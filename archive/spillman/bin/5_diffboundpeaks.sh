#!/usr/bin/env bash
#BSUB -J peakdiffs[1-2]
#BSUB -e peakdiffs.%J.%I.err
#BSUB -o peakdiffs.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P spillman

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/spillman/bin/config.sh

samples=(E21PR_plus E21PR_minus)
backgrounds=(E21PR_minus E21PR_plus)
sample=${samples[$(($LSB_JOBINDEX - 1))]}
background=${backgrounds[$(($LSB_JOBINDEX - 1))]}

getDifferentialPeaks $sample/peaks.txt $sample $background > $sample/peaks.diff.txt
