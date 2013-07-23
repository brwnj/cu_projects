#!/usr/bin/env bash
#BSUB -J findmotifs[1-8]
#BSUB -e findmotifs.%J.%I.err
#BSUB -o findmotifs.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 6
#BSUB -P thaikoottathil

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/thaikoottathil/bin/config.sh

tag=${TAGS[$(($LSB_JOBINDEX - 1))]}

findMotifsGenome.pl $RESULTS/$tag/peaks.txt hg19 $RESULTS/${tag}_motif -p 6 -size 50
