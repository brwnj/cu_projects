#! /usr/bin/env bash
#BSUB -J cyto_attrs
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

<<DOC
Gets the peak intensity for all miRNA seed hits to create a Cytoscape
attributes file for the edges.
DOC

set -o nounset -o pipefail -o errexit -x

# SAMPLES=(idx0 BT474herc MCF7 MCF7estr MDA231 BT474 BT474estr HS27A HS5 hMSC BMEC HUVEC)
# ALL_REPLICATES=(idx0 "PK23" "PK11 PK31 PK51" "PK12 PK32" "PK24 PK42 PK54" \
#             "PK21 PK41 PK52" "PK22 PK53" "MP1 MP21 MP35" "MP2 MP20 MP34"\
#             "MP36 MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA"\
#             "MP42.ACTG MP45.ACTG MP45.TCGA" "MP24 MP38")
# SAMPLE=${SAMPLES[$LSB_JOBINDEX]}
# REPLICATES=${ALL_REPLICATES[$LSB_JOBINDEX]}

SAMPLE=HELA
REPLICATES="helaa helab"

PROJECT=$HOME/projects/hits-clip
BIN=$PROJECT/bin
# GDARCHIVE=$PROJECT/data/combined.rmd.genomedata
GDARCHIVE=$HOME/projects/hits-clip/results/20120917/gd_20120917.rmd

CYTO_ATTR=$SAMPLE.intensity.eda
if [ ! -f $CYTO_ATTR ]; then
    echo "PeakIntensity" > $CYTO_ATTR

    STRANDS="pos neg"
    LENGTHS="7 8"
    for STRAND in $STRANDS; do
        for LENGTH in $LENGTHS; do
            TRACKS=$(for REP in $REPLICATES; do echo -n "$REP.$STRAND.rmd "; done)
            python $BIN/peak_intensity.py $SAMPLE.$STRAND.mirna.$LENGTH.notSIF.gz $GDARCHIVE $TRACKS >> $CYTO_ATTR
        done
    done
fi