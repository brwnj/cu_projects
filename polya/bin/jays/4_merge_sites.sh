#! /usr/bin/env bash
#BSUB -J merge_sites[1]
#BSUB -e merge_sites.%J.%I.err
#BSUB -o merge_sites.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

samples=(PC)
sample=${samples[$(($LSB_JOBINDEX - 1))]}
peaks="$RESULT/${sample}*/*.classified.bed.gz"

if [[ ! -d $RESULT/polya_sites ]]; then
    mkdir $RESULT/polya_sites
fi

sites13=$RESULT/polya_sites/$sample.sites.c13.bed.gz
sites1234=$RESULT/polya_sites/$sample.sites.c1234.bed.gz
siteswhole=$RESULT/polya_sites/$sample.sites.wholegene.bed.gz

if [[ ! -f $sites13 ]]; then
    python $BIN/merge_sites.py -n2 -c3 -c3a -x $XREF $EXONS $peaks \
        | bedtools sort -i - \
        | gzip -c > $sites13
fi
if [[ ! -f $sites1234 ]]; then
    python $BIN/merge_sites.py -n2 -c2 -c3 -c3a -c4 -x $XREF $EXONS $peaks \
        | bedtools sort -i - \
        | gzip -c > $sites1234
fi
if [[ ! -f $siteswhole ]]; then
    python $BIN/merge_sites.py -n2 -c2 -c3 -c3a -c4 -x $XREF $WHOLEGENE $peaks \
        | bedtools sort -i - \
        | gzip -c > $siteswhole
fi
