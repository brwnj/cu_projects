#! /usr/bin/env bash
#BSUB -J process_umi[1-5]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

<<DOC
remove duplicate umi entries at any given coordinate.
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(idx0 MP51 MP52 MP53 PK61 PK62)
SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}

BIN=$HOME/projects/hits-clip/bin

MOD=umi.unique

python $BIN/umitools-process-bam.py $SAMPLE.umi.bam $SAMPLE.$MOD.bam
samtools index $SAMPLE.$MOD.bam