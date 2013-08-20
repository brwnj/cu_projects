#!/usr/bin/env bash
set -o nounset

SAMPLES=(ACATCG GCCTAA TGGTCA)
NOVOIDX=$HOME/ref/sacCer3/sacCer3.nix
RESULTS=$HOME/projects/zhao/results/common
DATA=$HOME/projects/zhao/data/common
PICARD=$HOME/opt/picard
GATK=$HOME/opt/gatk/GenomeAnalysisTK.jar
REFERENCE=$HOME/ref/sacCer3/sacCer3.fa
