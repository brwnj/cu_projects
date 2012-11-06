#! /usr/bin/env bash
#BSUB -J index
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

DATA="$HOME/projects/zhao/data"
FASTA="$DATA/snRNA.fa"

bowtie2-build -f $FASTA snRNA