#!/usr/bin/env bash
#BSUB -J peakdiffs[1-20]
#BSUB -e peakdiffs.%J.%I.err
#BSUB -o peakdiffs.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

<<DOC
need to set this up differently...
DOC

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# if [[ $sample == *input* ]]; then
#     exit
# fi
#
# if [[ $sample == *mitotic* ]]; then
#     input=$RESULTS/mitotic_hela_input
# else
#     input=$RESULTS/hela_input
# fi
#
# peaks=$RESULTS/$sample/peaks.txt
# difpeaks=$RESULTS/$sample/differentially_bound_peaks.txt
# sampletagdir=$RESULTS/$sample
#
# getDifferentialPeaks $peaks $sampletagdir $input -F 2 > $difpeaks
