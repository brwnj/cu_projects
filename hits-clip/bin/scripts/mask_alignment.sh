#!/usr/bin/env bash
#BSUB -J novoalign[1-44]
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1

<<DOC
align masked fastqs using novoalign. toss anything that aligns in multiple
locations.
DOC

set -o nounset -o errexit -o pipefail -x

samples=(idx0 MP1 MP10 MP11 MP2 MP20
            MP21 MP22 MP23 MP24 MP30
            MP31 MP34 MP35 MP36 MP38
            MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG
            MP42.TCGA MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA
            MP45.ACTG MP45.TCGA MP7 MP9 PK11
            PK12 PK21 PK22 PK23 PK24
            PK31 PK32 PK33 PK41 PK42
            PK51 PK52 PK53 PK54)
sample=${samples[$LSB_JOBINDEX]}

masked_fastq=$HOME/projects/hits-clip/results/20130103/$sample.masked.fastq
novoidx=$HOME/projects/hits-clip/data/common/novoalign/hg18
align_result=$(dirname $masked_fastq)/$sample.masked.bam

novoalign -d $novoidx -f $masked_fastq -o SAM -r None \
    | samtools view -ShuF4 - \
    | samtools sort -o - $sample.temp -m 9500000000 > $align_result