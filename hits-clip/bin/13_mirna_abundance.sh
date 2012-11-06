#! /usr/bin/env bash
#BSUB -J mirna_abundance
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

<<DOC
From the miRNA coordinates, look up maximum peak intensity to represent the
abundance of each miRNA fragment.
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
# archive prior to removing the duplicates
# GDARCHIVE=$PROJECT/data/combined.genomedata
GDARCHIVE=$HOME/projects/hits-clip/results/20120917/gd_20120917
# bed was lifted over from hg19
REFERENCE=$PROJECT/data/common/mirbase/mirbase-18/hsa.hg18.bed.gz

USIF=$SAMPLE.network.unfiltered.sif
SIF=$SAMPLE.network.sif
CYTO_ATTR=$SAMPLE.abundance.noa
if [ ! -f $CYTO_ATTR ]; then
    # meeting cytoscape file/attribute requirements
    echo "mirnaAbundance" > $CYTO_ATTR

    python $BIN/mirna_abundance.py $REFERENCE $USIF $GDARCHIVE $REPLICATES >> $CYTO_ATTR
fi

if [ ! -f $SIF ]; then
    # filter the sif where abundance was found to be zero
    python $BIN/filter_sif.py $USIF $CYTO_ATTR > $SIF
fi