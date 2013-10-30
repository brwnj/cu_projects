#!/usr/bin/env bash
#BSUB -J align[1-14]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -q normal
#BSUB -n 8
#BSUB -P nicoli

<<DOC
align using novoalign
DOC

set -o nounset -o errexit -o pipefail -x
source $HOME/projects/nicoli/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq=$DATA/$sample.fastq.gz
results=$RESULTS/$sample
if [[ ! -d $results ]]; then
    mkdir $results
fi
bam=$results/$sample.bam
rmdups_bam=$results/$sample.rmd.bam

if [[ ! -f $bam ]]; then
    novoalign -d $NOVOIDX -f $fastq -a -o SAM -r A 20 -e 100 -s 2 -l 16 -c 8 -k \
        2> $sample.alignment_stats.txt \
        | samtools view -ShuF4 - \
        | samtools sort -o - $sample.temp -m 8G \
        > $bam
fi
if [[ ! -f $rmdups_bam ]]; then
    samtools rmdup -s $bam $rmdups_bam
fi
