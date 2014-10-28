#!/usr/bin/env bash
#BSUB -J "align[1-24]%8"
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 8
#BSUB -P leung

<<DOC
align rnaseq reads using gsnap
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/leung/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
r1=$DATA/common/${sample}_R1.fastq.gz
r2=$DATA/common/${sample}_R2.fastq.gz
results=$RESULTS/$sample
bam=$results/$sample.bam

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

if [[ ! -f $bam ]]; then
    gsnap -D $GMAPDB -d hg19_semiTotal --gunzip -t 8 -N 1 -s $SPLICESITES \
        --read-group-id $sample --read-group-name $sample --read-group-library $sample --read-group-platform Illumina \
        -B 3 -Q -n 20 -A sam $r1 $r2 \
        | samtools view -ShuF4 - \
        | samtools sort -o -m 8G - $sample.temp \
        > $bam
fi
