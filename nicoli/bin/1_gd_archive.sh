#!/usr/bin/env bash
#BSUB -J gdarch
#BSUB -e gdarch.%J.err
#BSUB -o gdarch.%J.out
#BSUB -q normal
#BSUB -P nicoli

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

bam2gd.py $SIZES $FASTAS $RESULTS/L*/*rmd.bam -o $GENOMEDATA -p nicoli
