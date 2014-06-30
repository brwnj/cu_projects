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

tempr1=${r1/.fastq.gz/.length_trimmed.fastq}
tempr2=${r2/.fastq.gz/.length_trimmed.fastq}

trimmedr1=${r1/.fastq.gz/.quality_trimmed.fastq.gz}
trimmedr2=${r2/.fastq.gz/.quality_trimmed.fastq.gz}

if [[ ! -f $trimmedr1 ]]; then
    trimfq.py -m 220 $r1 > $tempr1
    seqtk trimfq $tempr1 | gzip -c > $trimmedr1
    rm $tempr1
fi
if [[ ! -f $trimmedr2 ]]; then
    trimfq.py -m 220 $r2 > $tempr2
    seqtk trimfq $tempr2 | gzip -c > $trimmedr2
    rm $tempr1
fi
