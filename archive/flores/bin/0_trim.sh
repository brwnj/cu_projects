#! /usr/bin/env bash
#BSUB -J seqtk[1-4]
#BSUB -e seqtk.%J.%I.err
#BSUB -o seqtk.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P flores

<<DOC
quality trim the reads using seqtk
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/flores/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

r1=$DATA/${sample}_R1.fastq.gz
r2=$DATA/${sample}_R2.fastq.gz
trimmed_r1=$DATA/${sample}_R1.trm.fq.gz
trimmed_r2=$DATA/${sample}_R2.trm.fq.gz

seqtk trimfq $r1 | gzip -c > $trimmed_r1
seqtk trimfq $r2 | gzip -c > $trimmed_r2
