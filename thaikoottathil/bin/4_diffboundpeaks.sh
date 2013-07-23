#!/usr/bin/env bash
#BSUB -J peakdiffs[1-4]
#BSUB -e peakdiffs.%J.%I.err
#BSUB -o peakdiffs.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P thaikoottathil

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/thaikoottathil/bin/config.sh

samples=(STAT6_2h cJun_2h STAT6_24h cJun_24h)
backgrounds=(STAT6_2h_c cJun_2h_c STAT6_24h_c cJun_24h_c)
sample=$RESULTS/${samples[$(($LSB_JOBINDEX - 1))]}
background=$RESULTS/${backgrounds[$(($LSB_JOBINDEX - 1))]}
peaks=$sample/peaks.txt
diffpeaks=$sample/peaks.diff.txt

if [[ ! -f $diffpeaks ]] && [[ -f $peaks ]]; then
    getDifferentialPeaks $peaks $sample $background > $diffpeaks
fi
