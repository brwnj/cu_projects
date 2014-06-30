#! /usr/bin/env bash
#BSUB -J characterize_intronic_reads[1-44]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q short

<<DOC
Calculate stats on intronic reads from each sample
DOC

SAMPLEIDS=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
           MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG MP42.TCGA MP43.ACTG 
           MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11 
           PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53 
           PK54)
IDX=$(expr ${LSB_JOBINDEX} - 1)
SAMPLE=${SAMPLEIDS[$IDX]}

UTR3=$SAMPLE.intron.utr3.tab
UTR5=$SAMPLE.intron.utr5.tab
BIN=$HOME/projects/hits-clip/bin
REF=$HOME/projects/ref/hg18
DATA=$HOME/projects/hits-clip/results/common/samples
EXT=rmd.bed.gz

python $BIN/characterize_intronic_reads.py $DATA/$SAMPLE/$SAMPLE.$EXT \
        $REF/refseq.intron.bed.gz \
        $REF/refseq.utr3.bed.gz > $UTR3

python $BIN/characterize_intronic_reads.py $DATA/$SAMPLE/$SAMPLE.$EXT \
        $REF/refseq.intron.bed.gz \
        $REF/refseq.utr5.bed.gz > $UTR5