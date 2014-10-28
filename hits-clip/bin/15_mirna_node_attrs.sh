#! /usr/bin/env bash
#BSUB -J miRNA_node_features
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

<<DOC
Write miRNA's corresponding feature to node attributes file.
DOC

REF=$HOME/projects/ref/hg18
CDS=$REF/refseq.cds.bed.gz
EXON=$REF/refseq.exon.bed.gz
INTRON=$REF/refseq.intron.bed.gz
UTR3=$REF/refseq.utr3.bed.gz
UTR5=$REF/refseq.utr5.bed.gz

SRC=$HOME/projects/hits-clip/bin

MIRNA_BED=$HOME/projects/hits-clip/data/common/mirbase/mirbase-18/hsa.hg18.bed.gz

CYTO_ATTR=mirna.features.noa
echo "mirnaFeature" > $CYTO_ATTR

python $SRC/node_annotation.py --reference-beds \
        $CDS $EXON $INTRON $UTR3 $UTR5 \
        --names cds exon intron utr3 utr5 -- \
        $MIRNA_BED \
    | uniq \
    >> $CYTO_ATTR