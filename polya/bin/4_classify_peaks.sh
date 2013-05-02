#!/usr/bin/env bash
#BSUB -J classify_peaks[1-33]
#BSUB -e classify_peaks.%J.%I.err
#BSUB -o classify_peaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

<<DOC
classify peaks.
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/polya/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

result=$RESULT/$sample
peaks=$result/${sample}_peaks.bed.gz
posbg=$result/$sample.pos.bedgraph.gz
negbg=$result/$sample.neg.bedgraph.gz
classified=$result/${sample}_peaks.classified.bed.gz

if [[ ! -f $classified ]]; then
    pclassify.py $peaks $posbg $negbg $FASTA $CHROM_SIZES \
        | gzip -c > $classified
fi