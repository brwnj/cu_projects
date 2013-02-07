#!/usr/bin/env bash
#BSUB -J qvalues[1-24]
#BSUB -e qvalues.%J.%I.err
#BSUB -o qvalues.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1

<<DOC
get qvalues for all peaks
DOC

set -o nounset -o pipefail -o errexit -x

samples=(E1T1_Inf E1T1_Uninf E1T24_Inf E1T24_Uninf E1T2_Inf E1T2_Uninf
            E1T8_Inf E1T8_Uninf E2T1_Inf E2T1_Uninf E2T24_Inf E2T24_Uninf
            E2T2_Inf E2T2_Uninf E2T8_Inf E2T8_Uninf E3T1_Inf E3T1_Uninf
            E3T24_Inf E3T24_Uninf E3T2_Inf E3T2_Uninf E3T8_Inf E3T8_Uninf)
sample=${samples[$(expr $LSB_JOBINDEX - 1)]}

peaks=$HOME/projects/walter/results/common/$sample/$sample.novo.peaks.bed.gz
nullpeaks=$HOME/projects/walter/results/common/$sample/$sample.novo.peaks.shuffle.bed.gz
qvalues=$HOME/projects/walter/results/common/$sample/$sample.novo.qvalues.gz

peaktools-qvalues $realpeaks $nullpeaks | gzip -c > $qvalues