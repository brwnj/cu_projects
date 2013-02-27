#! /usr/bin/env bash
#BSUB -J novoalign[1-3]
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -q normal
#BSUB -n 8
#BSUB -P hits-clip

<<DOC
align using novoalign
DOC

set -o nounset -o errexit -o pipefail -x

# Samples for the job array
SAMPLEIDS=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
            MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG MP42.TCGA MP43.ACTG 
            MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11 
            PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53 
            PK54 helaa helab PK61 PK62 PK63)
SAMPLEIDS=(PK61 PK62 PK63)
SAMPLE=${SAMPLEIDS[$(($LSB_JOBINDEX - 1))]}

FASTQ=$HOME/projects/hits-clip/data/common/$(echo $SAMPLE | cut -d'.' -f1)/$SAMPLE.fastq.gz

RESULTS=$HOME/projects/hits-clip/results/common/samples/$SAMPLE
BAM=$RESULTS/$SAMPLE.bam
# duplicates removed
RBAM=$RESULTS/$SAMPLE.rmd.bam

NOVOIDX=$HOME/projects/hits-clip/data/common/novoalign/hg18

novoalign -d $NOVOIDX -f $FASTQ -a -o SAM -r A 20 -e 100 -s 2 -l 16 -c 8 \
    | samtools view -ShuF4 - \
    | samtools sort -o - $SAMPLE.temp -m 9500000000 \
    > $BAM

samtools rmdup -s $BAM $RBAM

# create bw for UCSC
SIZES=$HOME/ref/hg18/hg18.sizes
bam2bw $BAM $SIZES
bam2bw $RBAM $SIZES