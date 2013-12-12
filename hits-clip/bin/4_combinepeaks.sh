#!/usr/bin/env bash
#BSUB -J combinepeaks
#BSUB -e combinepeaks.%J.err
#BSUB -o combinepeaks.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_hitsclip

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/hits-clip/bin/config.sh

<<DOC
I run this inside the folder where 3_peaks.sh was run.
DOC

# serial processing of files... :(
for track in `genomedata-info tracknames_continuous $GENOMEDATA`; do
    sample=$(echo ${track  /_*} | cut -f1 -d.)
    peaks=$RESULTS/$sample/$track.peaks.bed
    shuffledpeaks=$RESULTS/$sample/$track.shuffle.peaks.bed
    rm -f $peaks $shuffledpeaks
    for (( i = 0; i < ${#CHROMOSOMES[@]}; i++ )); do
        chrom=${CHROMOSOMES[$i]}
        zcat $track.$chrom.peaks.bed.gz >> $peaks
        zcat $track.$chrom.shuffle.peaks.bed.gz >> $shuffledpeaks
    done
    bedSort $peaks $peaks
    bedSort $shuffledpeaks $shuffledpeaks
    gzip -f $peaks $shuffledpeaks
done
