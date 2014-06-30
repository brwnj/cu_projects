#!/usr/bin/env bash
#BSUB -J "novoalign[1-49]%5"
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -q normal
#BSUB -n 10
#BSUB -P hits-clip

<<DOC
align using novoalign
DOC

set -o nounset -o errexit -o pipefail -x
source $HOME/projects/hits-clip/bin/config.sh

sample=${SAMPLEIDS[$(($LSB_JOBINDEX - 1))]}
fastq=$DATA/$sample.fastq.gz
results=$RESULTS/$sample
bam=$results/$sample.bam
rmdups_bam=$results/$sample.rmd.bam
summary=$results/$sample.alignment_stats.txt

if [[ ! -d $results ]]; then
    mkdir $results
fi

if [[ ! -f $bam ]]; then
    novoalign -d $NOVOIDX -f $fastq -a -o SAM -r A 20 -e 100 -s 2 -l 16 -c 10 -k \
        2> $summary \
        | samtools view -ShuF4 - \
        | samtools sort -o - $sample.temp -m 8G \
        > $bam
fi

if [[ ! -f $rmdups_bam ]]; then
    samtools rmdup -s $bam $rmdups_bam
fi
