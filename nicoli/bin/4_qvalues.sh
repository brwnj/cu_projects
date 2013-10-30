#!/usr/bin/env bash
#BSUB -J qvalues[1-14]
#BSUB -e qvalues.%J.%I.err
#BSUB -o qvalues.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P nicoli

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
cutoff=.001

real=$RESULTS/$sample/$sample.rmd.peaks.bed.gz
null=$RESULTS/$sample/$sample.rmd.shuffle.peaks.bed.gz
allqvals=$RESULTS/$sample/$sample.rmd.peaks.qv.all.bed.gz
filteredqvals=$RESULTS/$sample/$sample.rmd.peaks.qv.passed_filter.bed.gz
summary=$RESULTS/$sample/qvalue_summary.txt

if [[ ! -f $filteredqvals ]]; then
    peaktools-qvalues -v $real $null | gzip -c > $allqvals
    awk -cheader -v CO=${cutoff} '$7<CO' $allqvals | gzip -c > $filteredqvals
    peaktools-qvalue-summary -t $cutoff $allqvals > $summary
    peaktools-qvalue-summary $allqvals >> $summary
fi
