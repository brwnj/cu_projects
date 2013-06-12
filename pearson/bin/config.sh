#!/usr/bin/env bash
set -o nounset -x

BASE=$HOME/projects/pearson
BIN=$BASE/bin
DATA=$BASE/data
RESULTS=$BASE/results
SAMPLES=(B1868_ATCACG_L001_R1_001 B1868_ATCACG_L001_R2_001
            DisA1_TGACCA_L001_R1_001 DisA1_TGACCA_L001_R2_001)
PAIRS=("$DATA/B1868_ATCACG_L001_R1_001.fastq.gz $DATA/B1868_ATCACG_L001_R2_001.fastq.gz"
        "$DATA/DisA1_TGACCA_L001_R1_001.fastq.gz $DATA/DisA1_TGACCA_L001_R2_001.fastq.gz")
NAMES=(B1868_ATCACG_L001 DisA1_TGACCA_L001)
# BOWTIE2_INDEX=$HOME/ref/tetrahymena/tta
# PICARD=$HOME/opt/picard-tools-1.79
REFERENCE=$HOME/ref/tetrahymena/Tetrahymena_thermophila.JCVI-TTA1-2.2.18.dna.toplevel.fa
GATK=$HOME/opt/gatk/GenomeAnalysisTK.jar
# GSNAPINDEX=$HOME/ref/gmapdb