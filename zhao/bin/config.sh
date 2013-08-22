#!/usr/bin/env bash
set -o nounset

# wt   ACATCG
# sup6 GCCTAA
# sup8 TGGTCA

SAMPLES=(ACATCG GCCTAA TGGTCA)
NOVOIDX=$HOME/ref/sacCer3/sacCer3.nix
RESULTS=$HOME/projects/zhao/results/common
DATA=$HOME/projects/zhao/data/common
PICARD=$HOME/opt/picard
GATK=$HOME/opt/gatk/GenomeAnalysisTK.jar
# need both for the reference fasta after sorting by chrom
# samtools faidx $REFERENCE
# java -jar ~/opt/picard/CreateSequenceDictionary.jar R=sacCer3.fa O=sacCer3.dict GENOME_ASSEMBLY=sacCer3 SPECIES=sacCer3
REFERENCE=$HOME/ref/sacCer3/sacCer3.fa
