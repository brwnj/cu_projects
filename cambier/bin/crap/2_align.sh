#! /usr/bin/env bash
#BSUB -J ohsnap.align[1-12]
#BSUB -e gsnap.%J.%I.err
#BSUB -o gsnap.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 4

<<DOC
align rnaseq samples using gsnap against mm9
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/cambier/bin/cambier.cfg

sample=${SAMPLES[$LSB_JOBINDEX]}
reads=$DATA/sample.fastq.gz
outdir=$RESULTS/common/$sample
bam=$outdir/$sample.bam
stats=$outdir/$sample.stats

if [ ! -f $bam ]; then
    gsnap -D $GMAPDB -d mm9 --gunzip -s $KNOWNSITES -v snp128_strict_wholeChrs \
        -n1 -B5 -Q --nofails -t4 -A sam --pairexpect=250 --pairdev=150 \
        --read-group-id=$sample --read-group-name=$sample $reads \
        | samtools view -ShuF 4 - \
        | samtools sort -o - $sample.temp -m 9500000000 > $bam
    samtools index $bam
fi
if [ ! -f $stats ]; then
    java -Xmx8g -jar $PICARD/CollectMultipleMetrics.jar \
        INPUT=$bam \
        REFERENCE_SEQUENCE=$REFERENCE \
        ASSUME_SORTED=true \
        OUTPUT=$stats \
        PROGRAM=CollectAlignmentSummaryMetrics \
        PROGRAM=CollectInsertSizeMetrics \
        PROGRAM=QualityScoreDistribution \
        PROGRAM=MeanQualityByCycle
fi