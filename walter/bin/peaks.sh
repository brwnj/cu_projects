#!/usr/bin/env bash
#BSUB -J peaktools[1-24]
#BSUB -e peaktools.%J.%I.err
#BSUB -o peaktools.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

<<DOC
call peaks using peaktools
DOC

set -o nounset -o pipefail -o errexit -x

samples=(idx0 E1T1_Inf E1T1_Uninf E1T24_Inf E1T24_Uninf E1T2_Inf E1T2_Uninf
            E1T8_Inf E1T8_Uninf E2T1_Inf E2T1_Uninf E2T24_Inf E2T24_Uninf
            E2T2_Inf E2T2_Uninf E2T8_Inf E2T8_Uninf E3T1_Inf E3T1_Uninf
            E3T24_Inf E3T24_Uninf E3T2_Inf E3T2_Uninf E3T8_Inf E3T8_Uninf)

sample=${samples[$LSB_JOBINDEX]}
gd=/vol1/home/brownj/projects/walter/results/common/genomedata_novo
peaks=/vol1/home/brownj/projects/walter/results/common/$sample/$sample.novo.peaks.bed
shuffledpeaks=/vol1/home/brownj/projects/walter/results/common/$sample/$sample.novo.peaks.shuffle.bed

strands=(pos neg)
symbols=(+ -)
for (( i = 0; i < 2; i++ )); do
    python ~/devel/peaktools/peaktools/identify_peaks.py -w 60 \
        -t $sample.novo_${strands[$i]} -s ${symbols[$i]} --trim-peaks $gd >> $peaks
    python ~/devel/peaktools/peaktools/identify_peaks.py -w 60 \
        -t $sample.novo_${strands[$i]} -s ${symbols[$i]} --shuffle-data $gd >> $shuffledpeaks
done
gzip $peaks $shuffledpeaks