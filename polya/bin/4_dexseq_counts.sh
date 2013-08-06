#! /usr/bin/env bash
#BSUB -J counts[1-63]
#BSUB -e counts.%J.%I.err
#BSUB -o counts.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sites=$RESULT/polya_sites/${sample:0:2}.sites.c13.bed.gz
results=$RESULT/$sample

for strand in pos neg; do
    bedg=$results/$sample.$strand.bedgraph.gz
    out=$results/$sample.$strand.counts.txt.gz
    # if [[ ! -f $out ]]; then
    python $BIN/read_counts.py $bedg $sites $CHROM_SIZES | gzip -c > $out
    # fi
done
