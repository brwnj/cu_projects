#!/usr/bin/env bash
#BSUB -J novoalign[1-24]
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 4

<<DOC
align reads that mapped to tb to hg19
DOC

set -o nounset -o pipefail -o errexit -x

samples=(idx0 E1T1_Inf E1T1_Uninf E1T24_Inf E1T24_Uninf E1T2_Inf E1T2_Uninf
            E1T8_Inf E1T8_Uninf E2T1_Inf E2T1_Uninf E2T24_Inf E2T24_Uninf
            E2T2_Inf E2T2_Uninf E2T8_Inf E2T8_Uninf E3T1_Inf E3T1_Uninf
            E3T24_Inf E3T24_Uninf E3T2_Inf E3T2_Uninf E3T8_Inf E3T8_Uninf)
sample=${samples[$LSB_JOBINDEX]}

novoidx=$HOME/ref/hg19/hg19.9606.novoidx
fastq=$HOME/projects/walter/data/20130204/$sample.tb.fq.gz
bam=$HOME/projects/walter/results/common/$sample/$sample.tb2hg19.novo.bam

novoalign -d $novoidx -f $fastq -o SAM -r A -s 2 -l 18 \
    -c 4 -k 2> $sample.tb2hg19.align_stats.txt \
    | samtools view -ShuF4 - \
    | samtools sort -o - $sample.temp -m 9500000000 > $bam