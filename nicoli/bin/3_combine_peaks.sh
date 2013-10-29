#!/usr/bin/env bash
#BSUB -J combinepeaks[1-10]
#BSUB -e combinepeaks.%J.%I.err
#BSUB -o combinepeaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P nicoli

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

peaks=$sample.rmd.peaks.bed
shuffledpeaks=$sample.rmd.shuffle.peaks.bed
rm -f $peaks $shuffledpeaks
for strand in pos neg; do
    for chrom in $CHROMOSOMES; do
        #regular peaks
        zcat $sample.$chrom.$strand.rmd.peaks.bed.gz >> $peaks
        #shuffled peaks
        zcat $sample.$chrom.$strand.rmd.shuffle.peaks.bed.gz >> $shuffledpeaks
    done
done
bedSort $peaks $peaks
bedSort $shuffledpeaks $shuffledpeaks

gzip -f $peaks $shuffledpeaks
peaks=$sample.rmd.peaks.bed.gz
shuffledpeaks=$sample.rmd.shuffle.peaks.bed.gz
cp $peaks $RESULTS/$sample/
cp $shuffledpeaks $RESULTS/$sample/
