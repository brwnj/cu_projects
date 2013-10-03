#!/usr/bin/env bash
#BSUB -J "novoalign[1-24]%5"
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 4

<<DOC
align to h37rv using novoalign
DOC

set -o nounset -o pipefail -o errexit -x

samples=(E1T1_Inf E1T1_Uninf E1T24_Inf E1T24_Uninf E1T2_Inf E1T2_Uninf
            E1T8_Inf E1T8_Uninf E2T1_Inf E2T1_Uninf E2T24_Inf E2T24_Uninf
            E2T2_Inf E2T2_Uninf E2T8_Inf E2T8_Uninf E3T1_Inf E3T1_Uninf
            E3T24_Inf E3T24_Uninf E3T2_Inf E3T2_Uninf E3T8_Inf E3T8_Uninf)
sample=${samples[$(($LSB_JOBINDEX - 1))]}

novoidx=$HOME/ref/tuberculosis/H37Rv.novoidx
# adapter are already trimmed
fastq=$HOME/projects/walter/data/20121005/$sample.trm.fq.gz
results=$HOME/projects/walter/results/common/$sample
bam=$results/$sample.tb.novo.bam

novoalign -d $novoidx -f $fastq -s 2 -l 17 -o SAM -r A \
    -c 4 -k 2> $sample.tb.align_stats.txt \
    | samtools view -ShuF4 - \
    | samtools sort -o - $sample.temp -m 9500000000 > $bam