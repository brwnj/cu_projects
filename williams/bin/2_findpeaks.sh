#!/usr/bin/env bash
#BSUB -J findpeaks[1-4]
#BSUB -e findpeaks.%J.%I.err
#BSUB -o findpeaks.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams/bin/config.sh

tag=${TAGS[$(($LSB_JOBINDEX - 1))]}
style=${STYLES[$(($LSB_JOBINDEX - 1))]}

input=$RESULTS/$tag
findPeaks $input -style $style -o auto -i $BACKGROUND
