#!/usr/bin/env bash
#BSUB -J maketagdir
#BSUB -e maketagdir.%J.%I.err
#BSUB -o maketagdir.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P cohrs_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/cohrs_chipseq/bin/config.sh

<<DOC
H3K9acB1_S5 H3K9ac  1
H3KacE1_S7  H3K9ac  2
H3PANA1_S1  H3PAN   1
H3PAND1_S3  H3PAN   2
INPUT_S2    input   1
INPUT_S4    input   2
RIGGC1_S6   R IgG   1
RIGGF1_S8   R IgG   2
DOC

echo "makeTagDirectory $RESULTS/H3K9ac $RESULTS/H3K9acB1_S5/H3K9acB1_S5.bam $RESULTS/H3KacE1_S7/H3KacE1_S7.bam" | bsez maketagdir -P cohrs_chipseq
echo "makeTagDirectory $RESULTS/H3PAN $RESULTS/H3PANA1_S1/H3PANA1_S1.bam $RESULTS/H3PAND1_S3/H3PAND1_S3.bam" | bsez maketagdir -P cohrs_chipseq
echo "makeTagDirectory $RESULTS/INPUT $RESULTS/INPUT_S2/INPUT_S2.bam $RESULTS/INPUT_S4/INPUT_S4.bam" | bsez maketagdir -P cohrs_chipseq
echo "makeTagDirectory $RESULTS/R_IgG $RESULTS/RIGGC1_S6/RIGGC1_S6.bam $RESULTS/RIGGC1_S8/RIGGC1_S8.bam" | bsez maketagdir -P cohrs_chipseq
