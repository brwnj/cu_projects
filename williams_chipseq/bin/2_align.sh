#!/usr/bin/env bash
#BSUB -J align[1-20]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>10] rusage[mem=10] span[hosts=1]"
#BSUB -n 8
#BSUB -P williams_chipseq

<<DOC
align reads using novoalign
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq=$DATA/$sample.trimmed.fastq.gz
results=$RESULTS/$sample

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

summary=$results/$sample.alignment_stats.txt
bam=$results/$sample.bam

if [[ ! -f $bam ]]; then
    novoalign -d $NOVOIDX -f $fastq -o SAM "@RG\tID:$sample\tSM:$sample\tPU:Illumina\tLB:PE" -r None -c 8 -k \
        2> $summary \
        | samtools view -ShuF4 - \
        | samtools sort -o -m 8G - $sample.temp \
        > $bam
fi
