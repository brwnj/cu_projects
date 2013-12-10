#! /usr/bin/env bash
#BSUB -J merge_sites[1-2]
#BSUB -e merge_sites.%J.%I.err
#BSUB -o merge_sites.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

samples=(MP PK)
sample=${samples[$(($LSB_JOBINDEX - 1))]}
peaks="$RESULTS/${sample}*/*.classified.bed.gz"

if [[ ! -d $POLYASITES ]]; then
    mkdir -p $POLYASITES
fi

testsites=$POLYASITES/$sample.test_sites.bed.gz
allsites=$POLYASITES/$sample.all_sites.bed.gz
siteswhole=$POLYASITES/$sample.wholegene_sites.bed.gz

slopsitesneg=$POLYASITES/$sample.test_sites.slop.$SLOP.neg.bed.gz
slopsitespos=$POLYASITES/$sample.test_sites.slop.$SLOP.pos.bed.gz

if [[ ! -d $POLYASITES ]]; then
    mkdir -p $POLYASITES
fi

if [[ ! -f $testsites ]]; then
    python $BIN/merge_sites.py -n2 -c3 -c3a -c5 -c5a -x $XREF $EXONS $peaks \
        | bedtools sort -i - \
        | gzip -c > $testsites
fi

if [[ ! -f $allsites ]]; then
    python $BIN/merge_sites.py -n2 -c2 -c3 -c3a -c4 -c5 -c5a -c6 -x $XREF $EXONS $peaks \
        | bedtools sort -i - \
        | gzip -c > $allsites
fi

if [[ ! -f $siteswhole ]]; then
    python $BIN/merge_sites.py -n2 -c2 -c3 -c3a -c4 -c5 -c5a -c6 -x $XREF $WHOLEGENE $peaks \
        | bedtools sort -i - \
        | gzip -c > $siteswhole
fi

# this adds slop to the polyA sites
if [[ ! -f $slopsitesneg ]]; then
    bedtools slop -b $SLOP -i $testsites -g $SIZES \
        | awk '$6 == "+"' \
        | bedtools sort -i - \
        | gzip -c > $slopsitesneg
fi
if [[ ! -f $slopsitespos ]]; then
    bedtools slop -b $SLOP -i $testsites -g $SIZES \
        | awk '$6 == "-"' \
        | bedtools sort -i - \
        | gzip -c > $slopsitespos
fi
