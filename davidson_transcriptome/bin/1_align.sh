#!/usr/bin/env bash
#BSUB -J align[1-10]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 8
#BSUB -P davidson_transcriptome

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/davidson_transcriptome/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq=$DATA/$sample.fastq.gz
results=$RESULTS/$sample
bam=$results/$sample.bam

if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $bam ]]; then
    novoalign -d $NOVOIDX -f $fastq -o SAM -r None -c 8 -k \
        2> $sample.alignment.txt \
        | samtools view -ShuF4 - \
        | samtools sort -o - $sample.temp -m 8G \
        > $bam
fi
