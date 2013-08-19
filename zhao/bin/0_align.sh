#!/usr/bin/env bash
#BSUB -J align[1-3]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P zhao

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/zhao/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
r1=$DATA/${sample}_R1.fastq.gz
r2=$DATA/${sample}_R2.fastq.gz
results=$RESULTS/$sample
bam=$results/$sample.bam
if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $bam ]]; then
    novoalign -d $NOVOIDX -f $r1 $r2 -o SAM "@RG\tID:$sample\tSM:$sample\tPU:Illumina\tLB:PE" -r None -i 250 100 -c 8 -k \
        2> $sample.alignment.txt \
        | samtools view -ShuF4 - \
        | samtools sort -o -m 8G - $sample.temp \
        > $bam
fi
