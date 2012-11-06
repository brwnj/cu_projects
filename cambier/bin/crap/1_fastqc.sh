#!/usr/bin/env bash
#BSUB -J fastqc[1-12]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

<<DOC
qc the reads
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/cambier/bin/cambier.cfg

sample=${SAMPLES[$LSB_JOBINDEX]}
outdir=$BASE/results/common/$sample
mkdir -p $outdir

$FASTQC --outdir $outdir $DATA/$sample.fastq.gz