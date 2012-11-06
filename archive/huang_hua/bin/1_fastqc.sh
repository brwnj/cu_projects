#!/usr/bin/env bash
#BSUB -J fastqc[1-6]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -n 4

<<DOC
qc the reads
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/huang_hua/bin/huang_hua.cfg

sample=${SAMPLES[${LSB_JOBINDEX}]}
outdir=$RESULTS/common/$sample/
mkdir -p $outdir

$FASTQC --outdir $outdir --threads 4 --format fastq $DATA/$sample.fastq.gz