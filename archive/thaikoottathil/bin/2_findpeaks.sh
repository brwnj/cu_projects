#!/usr/bin/env bash
#BSUB -J findpeaks[1-8]
#BSUB -e findpeaks.%J.%I.err
#BSUB -o findpeaks.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P thaikoottathil

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/thaikoottathil/bin/config.sh

tag=${TAGS[$(($LSB_JOBINDEX - 1))]}

input=$RESULTS/$tag
peaks=$input/peaks.txt
if [[ ! -f $peaks ]]; then
    findPeaks $input -style factor -o auto -i $BACKGROUND
fi
