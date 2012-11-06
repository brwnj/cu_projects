#! /usr/bin/env bash
#BSUB -J bowtie2.align
#BSUB -e bowtie2.%J.err
#BSUB -o bowtie2.%J.out
#BSUB -q normal
#BSUB -R "select[mem>30] rusage[mem=28] span[hosts=1]"
#BSUB -n 12

<<DOC
align fastqs using gsnap then sort and index the output.
DOC

set -o nounset -o pipefail -o errexit -x

DATA=/vol1/home/brownj/projects/peterson/data
READ1=$DATA/3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq
READ2=$DATA/3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq

GENOMEDIR=/vol1/home/gowank/Packages/GMAPDB

SAMPLE=3Poly_1794e3_Mutant
IDX=$HOME/projects/ref/mm9/bowtie2/mm9

# gsnap -D $GENOMEDIR -d mm9 -v snp128_strict_wholeChrs -n1 -B5 -Q \
#     --nofails -t12 -A sam --pairmax-dna=1000 $READ1 $READ2 > $SAM
# 
# samtools view -Shb $SAM > $UNSORTEDBAM
# samtools sort $UNSORTEDBAM $SAMPLE
# samtools index $BAM

# gatk requires read group to be set
bowtie2 -p12 -x $IDX -I 150 -X 1000 -q --very-sensitive --rg-id $SAMPLE \
    -1 $READ1 -2 $READ2 \
    | samtools view -bSh - > $SAMPLE.unsorted.bam
samtools sort $SAMPLE.unsorted.bam $SAMPLE
samtools index $SAMPLE.bam