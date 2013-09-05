#!/usr/bin/env bash
#BSUB -J classify_shifts[1-69]
#BSUB -e classify_shifts.%J.%I.err
#BSUB -o classify_shifts.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya
#BSUB -q short

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
infile=$DEXRESULTS/${sample}*
outfile=$SITESHIFTS/$sample.classified_shifts.txt.gz
sites=$RESULT/polya_sites/${sample:0:2}.sites.c13.bed.gz

# overwrite everything
python $BIN/classify_shifts.py $infile | gzip -c > $outfile

cd $SITESHIFTS
python $BIN/visualize_dexseq_shifts.py $outfile $sites
