#!/usr/bin/env bash
#BSUB -J findmotifs[1-4]
#BSUB -e findmotifs.%J.%I.err
#BSUB -o findmotifs.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 6
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams/bin/config.sh

tag=${TAGS[$(($LSB_JOBINDEX - 1))]}

peaks=$RESULTS/$tag/regions.txt
if [[ $tag == "AP2" ]]; then
    peaks=$RESULTS/$tag/peaks.txt
fi

findMotifsGenome.pl $peaks mm10 $RESULTS/${tag}_motif -p 6
