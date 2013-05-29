#! /usr/bin/env bash
#BSUB -J novoalign[1-4]
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P flores

<<DOC
align reads to <genome> using novoalign
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/flores/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

trimmed_r1=$DATA/${sample}_R1.trm.fq.gz
trimmed_r2=$DATA/${sample}_R2.trm.fq.gz

results=$RESULTS/$sample
stats=$results/$sample.alignment.hg19.txt
bam=$results/$sample.hg19.bam

if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $bam ]]; then
    novoalign -d $NOVOIDXHG19 -f $trimmed_r1 $trimmed_r2 -o SAM -r None -i 250 100 -c 8 -k \
        2> $stats \
        | samtools view -ShuF4 - \
        | samtools sort -o - $sample.temp -m 2G \
        > $bam
fi
