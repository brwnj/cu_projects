#!/usr/bin/env bash
#BSUB -J gdarch
#BSUB -e gdarch.%J.err
#BSUB -o gdarch.%J.out
#BSUB -q normal
#BSUB -P pillai_kabos_hitsclip

<<DOC
create genomedata archive.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/hits-clip/bin/config.sh

if [[ ! -d $GENOMEDATA ]]; then
    bam2gd.py $SIZES $FASTAS $bams -o $GENOMEDATA -p hitsclip &
fi
