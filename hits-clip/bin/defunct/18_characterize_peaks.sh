#! /usr/bin/env bash
#BSUB -J characterize_intronic_reads[1-10]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q short

<<DOC
Calculate stats on peaks of each case
DOC

SAMPLEIDS=(BT474herc MCF7 MCF7estr MDA231 BT474 BT474estr HS27A HS5 hMSC BMEC)
IDX=$(expr ${LSB_JOBINDEX} - 1)
SAMPLE=${SAMPLEIDS[$IDX]}

BIN=$HOME/projects/hits-clip/bin
REF=$HOME/projects/ref/hg18
DATA=$HOME/projects/hits-clip/results/common
EXT=peaks.trimmed.bed.gz
STATS=$SAMPLE.peaks.summary.txt

python $BIN/characterize_peaks.py $DATA/$SAMPLE/$SAMPLE.pos.$EXT \
        $DATA/$SAMPLE/$SAMPLE.neg.$EXT \
        $REF/refseq.intron.bed.gz intron -o $SAMPLE.intron > $STATS

python $BIN/characterize_peaks.py $DATA/$SAMPLE/$SAMPLE.pos.$EXT \
        $DATA/$SAMPLE/$SAMPLE.neg.$EXT \
        $REF/refseq.cds.bed.gz cds -o $SAMPLE.cds >> $STATS

python $BIN/characterize_peaks.py $DATA/$SAMPLE/$SAMPLE.pos.$EXT \
        $DATA/$SAMPLE/$SAMPLE.neg.$EXT \
        $REF/refseq.exon.bed.gz exon -o $SAMPLE.exon >> $STATS