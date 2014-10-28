#!/usr/bin/env bash
#BSUB -J findpeaks[1-20]
#BSUB -e findpeaks.%J.%I.err
#BSUB -o findpeaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

if [[ $sample == *input* ]]; then
    exit
fi

if [[ $sample == *mitotic* ]]; then
    input=$RESULTS/mitotic_hela_input
else
    input=$RESULTS/hela_input
fi

findPeaks $RESULTS/$sample -style factor -o auto -i $input
