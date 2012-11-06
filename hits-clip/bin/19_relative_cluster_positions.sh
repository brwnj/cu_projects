#! /usr/bin/env bash
#BSUB -J relative_cluster_positions[1-10]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q short

<<DOC
Output distribution of peak locations relative to feature
DOC

set -o nounset -o errexit -o pipefail -x

SAMPLEIDS=(BT474herc MCF7 MCF7estr MDA231 BT474 BT474estr HS27A HS5 hMSC BMEC)
IDX=$(expr ${LSB_JOBINDEX} - 1)
SAMPLE=${SAMPLEIDS[$IDX]}

BIN=$HOME/projects/hits-clip/bin
EXONS=$HOME/projects/ref/hg18/refseq.exons.fullname.bed
REF=$HOME/projects/ref/hg18/refseq.gtf.gz
DATA=$HOME/projects/hits-clip/results/common
EXT=peaks.trimmed.bed.gz

python $BIN/relative_cluster_positions.py $DATA/$SAMPLE/$SAMPLE.$EXT \
        $EXONS $REF stop_codon > $SAMPLE.peaklocs.stop_codon.txt
python $BIN/relative_cluster_positions.py $DATA/$SAMPLE/$SAMPLE.$EXT \
        $EXONS $REF start_codon > $SAMPLE.peaklocs.start_codon.txt