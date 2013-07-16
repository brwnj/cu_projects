#! /usr/bin/env bash
#BSUB -J callpeaks[1-3]
#BSUB -e callpeaks.%J.%I.err
#BSUB -o callpeaks.%J.%I.out
#BSUB -q normal
#BSUB -P hits-clip

# Samples for the job array $(ls ~brownj/projects/hits-clip/results/common/*)
# SAMPLEIDS=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
#             MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG MP42.TCGA MP43.ACTG
#             MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11
#             PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53
#             PK54 helaa helab PK61 PK62 PK63)
SAMPLEIDS=(PK61 PK62 PK63)
SAMPLE=${SAMPLEIDS[$(($LSB_JOBINDEX - 1))]}
#PK61.neg.rmd
GENOMEDATA=$HOME/projects/hits-clip/data/combined.rmd.genomedata
STRANDS="pos neg"
SYMBOLS=(+ -)
COUNT=0

for STRAND in $STRANDS; do
    # regular peaks
    RUNSCRIPT="$SAMPLE.$STRAND.peaks.sh"
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J $SAMPLE.$STRAND.peaks[1-25]" >> $RUNSCRIPT
    echo "#BSUB -e $SAMPLE.$STRAND.%J.%I.err" >> $RUNSCRIPT
    echo "#BSUB -o $SAMPLE.$STRAND.%J.%I.out" >> $RUNSCRIPT
    echo "#BSUB -q normal" >> $RUNSCRIPT
    echo "#BSUB -P hits-clip" >> $RUNSCRIPT
    echo "
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
CHROMOSOME=\${CHROMOSOMES[\$((\$LSB_JOBINDEX - 1))]}
peaktools-identify-peaks -t $SAMPLE.$STRAND.rmd -c \$CHROMOSOME -w 30 -v -s ${SYMBOLS[$COUNT]} $GENOMEDATA | gzip -c > $SAMPLE.\$CHROMOSOME.$STRAND.rmd.peaks.bed.gz
" >> $RUNSCRIPT
    bsub < $RUNSCRIPT

    #shuffled peaks
    RUNSCRIPT="$SAMPLE.$STRAND.shuffle.peaks.sh"
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J $SAMPLE.$STRAND.shuffle[1-25]" >> $RUNSCRIPT
    echo "#BSUB -e %J.%I.err" >> $RUNSCRIPT
    echo "#BSUB -o %J.%I.out" >> $RUNSCRIPT
    echo "#BSUB -q normal" >> $RUNSCRIPT
    echo "#BSUB -P hits-clip" >> $RUNSCRIPT
    echo "
CHROMOSOMES=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
CHROMOSOME=\${CHROMOSOMES[\$((\$LSB_JOBINDEX - 1))]}
peaktools-identify-peaks -t $SAMPLE.$STRAND.rmd -c \$CHROMOSOME -w 30 -v -s ${SYMBOLS[$COUNT]} --shuffle-data $GENOMEDATA | gzip -c > $SAMPLE.\$CHROMOSOME.$STRAND.rmd.shuffle.peaks.bed.gz
" >> $RUNSCRIPT
    bsub < $RUNSCRIPT

    COUNT=$(expr $COUNT + 1)
done