#!/usr/bin/env bash
#BSUB -J read_stats
#BSUB -e read_stats.%J.err
#BSUB -o read_stats.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

ONE=""
TWO=""
THREE=""
FOUR=""
STATS=$READS/count_table.txt

# have some extra barcodes and only want to compile stats on those specified
# in SAMPLES
for f in ${SAMPLES[@]}; do
    ONE=$ONE" "$READS/${f}_R1.fastq.gz
    TWO=$TWO" "$READS/${f}_R1.umifiltered.fastq.gz
    THREE=$THREE" "$READS/${f}_R1.trimmed.fastq.gz
    FOUR=$FOUR" "$READS/$f.joined.fastq.gz
done

python ~/projects/bennett/bin/read_counts.py --one $ONE --two $TWO --three $THREE --four $FOUR > $STATS
