#! /usr/bin/env bash
#BSUB -J combinepeaks[1-46]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

SAMPLEIDS=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
           MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG MP42.TCGA MP43.ACTG 
           MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11 
           PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53 
           PK54 helaa helab)
IDX=$(expr $LSB_JOBINDEX - 1)
SAMPLE=${SAMPLEIDS[$IDX]}

STRANDS="pos neg"
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM"
DATA=/vol1/home/brownj/projects/hits-clip/results/common/samples
COMBINEDBED=$SAMPLE.peaks.bed
SHUFFLECOMBINED=$SAMPLE.shuffle.peaks.bed
for STRAND in $STRANDS; do
    for CHROM in $CHROMOSOMES; do
        #regular peaks
        zcat $SAMPLE.$CHROM.$STRAND.peaks.bed.gz >> $COMBINEDBED
        #shuffled peaks
        zcat $SAMPLE.$CHROM.$STRAND.shuffle.peaks.bed.gz >> $SHUFFLECOMBINED
    done
done
bedSort $COMBINEDBED $COMBINEDBED
bedSort $SHUFFLECOMBINED $SHUFFLECOMBINED
gzip -f $COMBINEDBED $SHUFFLECOMBINED
cp $COMBINEDBED.gz $DATA/$SAMPLE/$COMBINEDBED.gz
cp $SHUFFLECOMBINED.gz $DATA/$SAMPLE/$SHUFFLECOMBINED.gz