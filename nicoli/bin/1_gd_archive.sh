#!/usr/bin/env bash
#BSUB -J gdarch
#BSUB -e gdarch.%J.err
#BSUB -o gdarch.%J.out
#BSUB -q normal
#BSUB -P nicoli

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

bams=$RESULTS/L*/*bam
for file in $bams; do
    samtools index $file
done

bam2gd.py $SIZES $FASTAS $bams -o $GENOMEDATA -p nicoli
