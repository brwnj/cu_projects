#!/usr/bin/env bash
#BSUB -J findmotifs[1-2]
#BSUB -e findmotifs.%J.%I.err
#BSUB -o findmotifs.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 12
#BSUB -P spillman

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/spillman/bin/config.sh

samples=(E21PR_plus E21PR_minus)
sample=${samples[$(($LSB_JOBINDEX - 1))]}

findMotifsGenome.pl $sample/peaks.txt hg19 ${sample}_motif -p 12 -size 50
