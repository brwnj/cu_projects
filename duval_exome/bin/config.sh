#!/usr/bin/env bash
SAMPLES=(1025 113 166408 172160 22 353GAAACC 353GGCTAC 36 400 500 522 730 735 868 Bililey KTCC)
RESULTS=$HOME/projects/duval_exome/results/common
DATA=$HOME/projects/duval_exome/data/common
NOVOIDX=$HOME/ref/canFam3/canFam3.novoidx
PICARD=$HOME/opt/picard
GATK=$HOME/opt/gatk/GenomeAnalysisTK.jar
REFERENCE=$HOME/ref/canFam3/canFam3.fa
SNPEFF=$HOME/opt/snpeff/snpEff.jar
SNPSIFT=$HOME/opt/snpeff/SnpSift.jar
SNPEFFCONFIG=$HOME/opt/snpeff/snpEff.config
BIN=$HOME/projects/duval_exome/bin
ANNOTATION=$HOME/ref/canFam3/canFam3.descriptions.gz
TARGETS=$HOME/ref/canFam3/canFam3.SureDesign.bed.gz
NORMALVCF=$RESULTS/normal.vcf
HUB=$RESULTS/hub
SIZES=$HOME/ref/canFam3/canFam3.sizes
