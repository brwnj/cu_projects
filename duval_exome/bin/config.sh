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

declare -A HUBCOLORS
HUBCOLORS=(
[1025]="35,139,69"
[113]="136,65,157"
[166408]="43,140,190"
[172160]="215,48,31"
[22]="5,11,176"
[353GAAACC]="2,129,138"
[353GGCTAC]="206,18,86"
[36]="174,1,126"
[400]="35,132,67"
[500]="34,94,168"
[522]="204,76,2"
[730]="227,26,28"
[735]="102,194,164"
[868]="140,150,198"
[Bililey]="123,204,196"
[KTCC]="252,141,89"
)
