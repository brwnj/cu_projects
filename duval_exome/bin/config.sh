#!/usr/bin/env bash
SAMPLES=(1025 113 166408 172160 22 353GAAACC 353GGCTAC 36 400 500 522 730 735 868 Bililey KTCC)
RESULTS=$HOME/projects/duval_exome/results/common
DATA=$HOME/projects/duval_exome/data/common
NOVOIDX=$HOME/ref/canis_familiaris/canis_familiaris.masked.nix
PICARD=$HOME/opt/picard
GATK=$HOME/opt/gatk/GenomeAnalysisTK.jar
REFERENCE=$HOME/ref/canis_familiaris/canis_familiaris.fa
SNPEFF=$HOME/opt/snpeff/snpEff.jar
SNPEFFCONFIG=$HOME/opt/snpeff/snpEff.config
BIN=$HOME/projects/duval_exome/bin
ANNOTATION=$HOME/ref/canis_familiaris/canfam3.descriptions.gz
