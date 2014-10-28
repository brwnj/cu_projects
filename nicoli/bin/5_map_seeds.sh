#!/usr/bin/env bash
#BSUB -J seedmapping[1-2]
#BSUB -e seedmapping.%J.%I.err
#BSUB -o seedmapping.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P nicoli

<<DOC

DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
results=$RESULTS/$sample
pospeak=$results/$sample.rmd_pos.peaks.qv.passed_filter.bed.gz
negpeak=$results/$sample.rmd_neg.peaks.qv.passed_filter.bed.gz

posout=$results/$sample.rmd_pos.seed7_mapping.txt.gz
negout=$results/$sample.rmd_neg.seed7_mapping.txt.gz

if [[ ! -f $SEED7 ]]; then
    cliptools-report-mirna-seeds -s 2 -e 8 $MIRBASE | gzip -c > $SEED7
fi
if [[ ! -f $SEED8 ]]; then
    cliptools-report-mirna-seeds -s 2 -e 9 $MIRBASE | gzip -c > $SEED8
fi

echo "cliptools-aggregate-peaks $pospeak $GENEANNOTATION $SEED7 $GENOMEDATA | gzip -c > $posout" | bsez seed2peak -P $PI
echo "cliptools-aggregate-peaks $negpeak $GENEANNOTATION $SEED7 $GENOMEDATA | gzip -c > $negout" | bsez seed2peak -P $PI

posout=$results/$sample.rmd_pos.seed8_mapping.txt.gz
negout=$results/$sample.rmd_neg.seed8_mapping.txt.gz

echo "cliptools-aggregate-peaks $pospeak $GENEANNOTATION $SEED8 $GENOMEDATA | gzip -c > $posout" | bsez seed2peak -P $PI
echo "cliptools-aggregate-peaks $negpeak $GENEANNOTATION $SEED8 $GENOMEDATA | gzip -c > $negout" | bsez seed2peak -P $PI
