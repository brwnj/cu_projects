#!/usr/bin/env bash
#BSUB -J alignment[1-2]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>21] rusage[mem=20] span[hosts=1]"
#BSUB -n 4

<<DOC
align paired-end samples using gsnap
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/pearson/bin/base.cfg

pair=${PAIRS[$LSB_JOBINDEX]}
name=${NAMES[$LSB_JOBINDEX]}

outdir=$RESULTS/common/$name
bam=$outdir/$name.mic.bam
stats=$outdir/$name.mic.stats
# mkdir -p $outdir

if [ ! -f $bam ]; then
    gsnap --gunzip -d tta_mic -D $GSNAPINDEX --input-buffer-size=100000 -B5 -t4 -A sam --nofails \
        -n1 --pairexpect=250 --pairdev=150 --read-group-id=$name \
        --read-group-name=$name $pair \
        | samtools view -ShuF 4 - | samtools sort -o - $name.temp -m 2G > $bam
    samtools index $bam
fi
# if [ ! -f $stats ]; then
#     java -Xmx8g -jar $PICARD/CollectMultipleMetrics.jar \
#         INPUT=$bam \
#         REFERENCE_SEQUENCE=$REFERENCE \
#         ASSUME_SORTED=true \
#         OUTPUT=$stats \
#         PROGRAM=CollectAlignmentSummaryMetrics \
#         PROGRAM=CollectInsertSizeMetrics \
#         PROGRAM=QualityScoreDistribution \
#         PROGRAM=MeanQualityByCycle
# fi