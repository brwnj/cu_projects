#! /usr/bin/env bash
#BSUB -J bg.bw[1-69]
#BSUB -e bg.bw.%J.%I.err
#BSUB -o bg.bw.%J.%I.out
#BSUB -q short
#BSUB -P pillai_kabos_polya

<<DOC
Convert aligned BAMs to bedgraph and bigwig format. Overwrites any existing files.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULTS/$sample
umibam=$results/$sample.UMIs_not_removed.bam
bam=$results/$sample.bam

bam2bw.py -5 -b -v $bam $SIZES pillai_kabos_polya
bam2bw.py -5 -b -v $umibam $SIZES pillai_kabos_polya
