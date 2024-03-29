#!/usr/bin/env bash
#BSUB -J novoalign[1-6]
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 4

<<DOC
realign only to mirbase
DOC

set -o nounset -o pipefail -o errexit -x

samples=(idx0 MDX_22 MDX_23 MDX_24 WT_21 WT_25 WT_42)
sample=${samples[$LSB_JOBINDEX]}

novoidx=$HOME/ref/mirbase/19/mature.mmu.novoidx
fastq=$HOME/projects/leinwand/data/20121101/$sample.trim.fastq.gz
results=$HOME/projects/leinwand/results/common/$sample
bam=$results/$sample.bam

novoalign -d $novoidx -f $fastq -s 2 -l 16 -o SAM -r A 50 \
    -c 4 -k 2> $sample.align_stats.txt \
    | samtools view -ShuF4 - \
    | samtools sort -o - $sample.temp -m 9500000000 > $bam