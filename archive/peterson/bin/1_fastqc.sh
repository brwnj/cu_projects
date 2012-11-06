#! /usr/bin/env bash
#BSUB -J fastqc
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

<<DOC
qc the reads
DOC

PROJECT=$HOME/projects/peterson

$HOME/opt/FastQC/fastqc --outdir=$(pwd -P) \
    $PROJECT/data/3Poly_1794e3_Mutant_NoIndex_L001_R1_001.fastq \
    $PROJECT/data/3Poly_1794e3_Mutant_NoIndex_L001_R2_001.fastq