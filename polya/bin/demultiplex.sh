#!/usr/bin/env bash
#BSUB -J demultiplex[1-2]
#BSUB -e demultiplex.%J.%I.err
#BSUB -o demultiplex.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 6

<<DOC
demultiplex using ea-tools' fastq-multx
DOC

set -o nounset -o errexit -x

SAMPLES=(idx0 6 7)
SAMPLE=${SAMPLES[$LSB_JOBINDEX]}

barcodes=$HOME/projects/polya/data/20130104/L00${SAMPLE}.txt
fastq=$HOME/projects/polya/results/20130110/${SAMPLE}.fq.gz

fastq-multx -B $barcodes -m 1 $fastq -o %.fq