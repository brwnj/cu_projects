#! /usr/bin/env bash
#BSUB -J celera[1-7]
#BSUB -e celera.%J.%I.err
#BSUB -o celera.%J.%I.out
#BSUB -q bigmem
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -P santangelo

<<DOC
Run Celera on all samples.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/santangelo/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
frg=$READS/$sample.frg
out=$RESULTS/$sample

asm=$out/$sample.asm

if [[ ! -f $asm ]]; then
    runCA -d $out -p $sample $frg
fi