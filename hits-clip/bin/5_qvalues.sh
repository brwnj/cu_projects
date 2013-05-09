#!/usr/bin/env bash
#BSUB -J qvalues[1-3]
#BSUB -e qvalues.%J.%I.err
#BSUB -o qvalues.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_hitsclip

set -o nounset -o pipefail -o errexit -x

# SAMPLEIDS=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
#            MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG MP42.TCGA MP43.ACTG 
#            MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11 
#            PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53 
#            PK54)
SAMPLEIDS=(PK61 PK62 PK63)
SAMPLE=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
DATA=$HOME/projects/hits-clip/results/common/samples
CUTOFF=.001

REAL=$DATA/$SAMPLE/$SAMPLE.peaks.rmd.bed.gz
NULL=$DATA/$SAMPLE/$SAMPLE.shuffle.peaks.rmd.bed.gz
QVALS=$DATA/$SAMPLE/$SAMPLE.peaks.rmd.qv.bed.gz
VALIDPEAKS=$DATA/$SAMPLE/$SAMPLE.peaks.rmd.qv$CUTOFF.bed.gz
SUMMARY=$DATA/$SAMPLE/qvalue_summary.txt

peaktools-qvalues -v $REAL $NULL | gzip -c > $QVALS
bioawk -c header -v CO=${CUTOFF} '$7<CO' $QVALS | gzip -c > $VALIDPEAKS
peaktools-qvalue-summary -t $CUTOFF $QVALS > $SUMMARY
peaktools-qvalue-summary $QVALS >> $SUMMARY
