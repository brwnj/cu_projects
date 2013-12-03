#!/usr/bin/env bash
#BSUB -J preprocess_fq[1-8]
#BSUB -e preprocess_fq.%J.%I.err
#BSUB -o preprocess_fq.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P cohrs_chipseq

<<DOC
trim the fastq using seqtk
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/cohrs_chipseq/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1=$DATA/${sample}_R1.fastq.gz
r2=$DATA/${sample}_R2.fastq.gz

tempr1=
tempr2=

trimmedr1=
trimmedr2=

python trimfq -m 220

seqtk trimfq -e 80 $r1 > 
seqtk trimfq -e 80 $r2
