#! /usr/bin/env bash
#BSUB -J identify_null_peaks[1-6]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

set -o nounset -o errexit -x

samples=(idx0
MDX_22_AGTTCC_L003_R1_001
MDX_23_ATGTCA_L003_R1_001
MDX_24_CCGTCC_L003_R1_001
WT_21_AGTCAA_L003_R1_001
WT_25_GTAGAG_L003_R1_001
WT_42_GTCCGC_L003_R1_001)

sample=${samples[${LSB_JOBINDEX}]}
strands="pos neg"
symbols="+ -"

results=/vol1/home/brownj/projects/leinwand/results/common
genomedata=$results/genomedata

for strand in $strands; do
    for symbol in $symbols; do
        peaktools-identify-peaks -v -t ${sample}_${strand} -w 50 -s $symbol \
        --trim-peaks --shuffle-data $genomedata \
        > $results/$sample/$sample.$strand.null.peaks
    done
done