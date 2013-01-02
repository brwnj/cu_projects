#! /usr/bin/env bash
<<DOC
create genomedata archive with sequence data
DOC

set -o nounset -o pipefail -o errexit -x

CHROM=chr1
SAMPLEIDS="E1T1_Inf E1T1_Uninf E1T24_Inf E1T24_Uninf E1T2_Inf E1T2_Uninf
            E1T8_Inf E1T8_Uninf E2T1_Inf E2T1_Uninf E2T24_Inf E2T24_Uninf
            E2T2_Inf E2T2_Uninf E2T8_Inf E2T8_Uninf E3T1_Inf E3T1_Uninf
            E3T24_Inf E3T24_Uninf E3T2_Inf E3T2_Uninf E3T8_Inf E3T8_Uninf"
STRANDS="pos neg"
DATA=/vol1/home/brownj/projects/walter/results/common

VARS=vars.$CHROM.sh
echo "SEQDIR=/vol1/home/brownj/ref/tuberculosis" > $VARS
echo "GENOMEDATADIR=$CHROM.gd" >> $VARS
echo "BED_TRACKNAMES=($(for SAMPLE in $SAMPLEIDS; do for STRAND in $STRANDS; do echo -n "$SAMPLE.$STRAND "; done; done))" >> $VARS
echo "BED_FILENAMES=($(for SAMPLE in $SAMPLEIDS; do for STRAND in $STRANDS; do echo -n "$DATA/$SAMPLE/$SAMPLE.H37Rv.$STRAND.bedgraph.gz "; done; done))" >> $VARS

JOBNAME="genomedata_load.$CHROM"
bsub -J $JOBNAME -o %J.out -e %J.err "load_chr_sequence.sh $VARS $CHROM"