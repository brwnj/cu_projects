#!/usr/bin/env bash
#BSUB -J findmotifs[1-20]
#BSUB -e findmotifs.%J.%I.err
#BSUB -o findmotifs.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 10
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
peaks=$RESULTS/$sample/peaks.txt
out=$RESULTS/$sample

if [[ $sample == *input* ]]; then
    exit
fi

findMotifsGenome.pl $peaks hg19 $out -p 10 -size 200
