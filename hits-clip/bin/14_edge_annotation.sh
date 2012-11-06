#! /usr/bin/env bash
#BSUB -J mirna_feature
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

<<DOC
From the miRNA coordinates (taken from .notSIF), annotate the feature of the 
edge.
DOC

set -o nounset -o errexit -x

# SAMPLES=(idx0 BT474herc MCF7 MCF7estr MDA231 BT474 BT474estr HS27A HS5 hMSC BMEC HUVEC)
# ALL_REPLICATES=(idx0 "PK23" "PK11 PK31 PK51" "PK12 PK32" "PK24 PK42 PK54" \
#             "PK21 PK41 PK52" "PK22 PK53" "MP1 MP21 MP35" "MP2 MP20 MP34"\
#             "MP36 MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA"\
#             "MP42.ACTG MP45.ACTG MP45.TCGA" "MP24 MP38")
# SAMPLE=${SAMPLES[$LSB_JOBINDEX]}
# REPLICATES=${ALL_REPLICATES[$LSB_JOBINDEX]}

SAMPLE=HELA
REPLICATES="helaa helab"

PROJECT=/vol1/home/brownj/projects/hits-clip
SRC=$PROJECT/bin

# archive prior to removing the duplicates
# GDARCHIVE=$PROJECT/data/combined.genomedata
GDARCHIVE=$HOME/projects/hits-clip/results/20120917/gd_20120917

REF=$HOME/ref/hg18
CDS=$REF/refseq.cds.bed.gz
EXON=$REF/refseq.exon.bed.gz
INTRON=$REF/refseq.intron.bed.gz
UTR3=$REF/refseq.utr3.bed.gz
UTR5=$REF/refseq.utr5.bed.gz

CYTO_ATTR=$SAMPLE.feature.eda
if [ ! -f $CYTO_ATTR ]; then
    echo "peakFeature" > $CYTO_ATTR

    STRANDS="pos neg"
    LENGTHS="7 8"
    for STRAND in $STRANDS; do
        for LENGTH in $LENGTHS; do
            python $SRC/peak_annotation.py \
                --beds $CDS $EXON $INTRON $UTR3 $UTR5 \
                --names cds exon intron utr3 utr5 -- \
                $SAMPLE.$STRAND.mirna.$LENGTH.notSIF.gz \
            >> $CYTO_ATTR
        done
    done
fi