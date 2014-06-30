#!/usr/bin/env bash
#BSUB -J qc[1-6]
#BSUB -e qc.%J.%I.err
#BSUB -o qc.%J.%I.out
#BSUB -n 4
#BSUB -q normal

<<DOC
align the reads using rum
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/davidson/bin/davidson.cfg

sample=${SAMPLES[$LSB_JOBINDEX]}
shortname=$(echo $sample | cut -f1 -d_)
outdir=$RESULTS/common/$shortname
date=20120924
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

read1=$(echo $sample | cut -f1 -d" ")
read2=$(echo $sample | cut -f2 -d" ")

# quality control the reads
$FASTQC --outdir $outdir --threads 4 --format fastq $DATA/$date/$read1 $DATA/$date/$read2