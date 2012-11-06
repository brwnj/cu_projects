#! /usr/bin/env bash
#BSUB -J index
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -R "rusage[swp=2] select[swp>2]"
#BSUB -q normal

FASTA=/vol1/home/brownj/projects/zhao/data/snRNA.fa

novoindex snRNA $FASTA