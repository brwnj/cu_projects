#!/usr/bin/env bash
set -o nounset -x

BASE=$HOME/projects/pearson
BIN=$BASE/bin
DATA=$BASE/data
RESULTS=$BASE/results
FASTQC=$HOME/opt/fastqc/fastqc
SAMPLES=(idx0
B1868_ATCACG_L001_R1_001
B1868_ATCACG_L001_R2_001
DisA1_TGACCA_L001_R1_001
DisA1_TGACCA_L001_R2_001)
PAIRS=(idx0
"$DATA/B1868_ATCACG_L001_R1_001.fastq.gz $DATA/B1868_ATCACG_L001_R2_001.fastq.gz"
"$DATA/DisA1_TGACCA_L001_R1_001.fastq.gz $DATA/DisA1_TGACCA_L001_R2_001.fastq.gz")
NAMES=(idx0
B1868_ATCACG_L001
DisA1_TGACCA_L001)
BOWTIE2_INDEX=$HOME/ref/tetrahymena/tta
PICARD=$HOME/opt/picard-tools-1.79
# REFERENCE=$HOME/ref/tetrahymena/tta1_oct2008_finalrelease.assembly.fa
REFERENCE=$HOME/ref/tetrahymena/tetrahymena_thermophila_sb210__mic__2_supercontigs.fasta
GATK=$HOME/opt/gatk/GenomeAnalysisTK.jar
GSNAPINDEX=$HOME/ref/gmapdb