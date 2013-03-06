#! /usr/bin/env bash
#BSUB -J novoalign[1-4]
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P duval_chipseq

<<DOC
Align using Novoalign, suppressing all reads that align more than once.
DOC

set -o nounset -o pipefail -o errexit -x

samples=(1 2 3 4)
# yes, i know...
sample=${samples[$(($LSB_JOBINDEX - 1))]}

fastq=$HOME/projects/duval/data/20130215/$sample.trm.fq.gz
novoidx=$HOME/ref/hg19/hg19.9606.novoidx
results=$HOME/projects/duval/results/common/$sample
bam=$results/$sample.bam

if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $bam ]]; then
novoalign -d $novoidx -f $fastq -o SAM -r None -l 20 -s 3 -c 8 -k \
    2> $sample.align_stats.txt \
    | samtools view -ShuF4 - \
    | samtools sort -o - $sample.temp -m 9500000000 \
    > $bam
fi
samtools index $bam