#!/usr/bin/env bash
set -o nounset

SAMPLES=(
130-4_20h
130-4-HSV-1_20h
143_1_20h
143_1HSV1_20h
152_6_20h
152_6hsv1_20h
169-2_20h
169-2HSV-1_20h
17-5HSV-1_24h
17-5sham_24h
24-3-HSV-1_24h
24-3sham_24h
25_9HSV1_24h
25_9sham_24h
329-5_20h
329-5hsv-1_20h
35-5HSV-1_21h
35-5sham_21h
373-7Rhsv-1_21h
373-7Rsham_21h
58-1_21h
58-1hsv-1_21h
88-3_24h
88-3hsv-1_24h
)
RAWDATA=$HOME/projects/leung/data/20131127
DATA=$HOME/projects/leung/data
RESULTS=$HOME/projects/leung/results/common
SPLICESITES=~analysiscore/genomes/GMAPDB/hg19_semiTotal/hg19_semiTotal.maps/hg19_semiTotal.ensembl.splicesites
GMAPDB=~analysiscore/genomes/GMAPDB/hg19_semiTotal
REFERENCE=$HOME/ref/hg19/hg19.fa
PICARD=$HOME/opt/picard
SNPEFF=$HOME/opt/snpeff/snpEff.jar
SNPSIFT=$HOME/opt/snpeff/SnpSift.jar
SNPEFFCONFIG=$HOME/opt/snpeff/snpEff.config
HUB=$RESULTS/hub
