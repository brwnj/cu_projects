#! /usr/bin/env bash
#BSUB -J bowtie[1-5]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

<<DOC
Align using Bowtie, suppressing all reads that align more than once.
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(idx0 MP51 MP52 MP53 PK61 PK62)
SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}

BIDX=/vol3/home/jhessel/projects/bowtie/indices/hg18

READS=/vol1/home/brownj/projects/polya/data/20121008

gunzip $READS/$SAMPLE.umi.fq.gz
# bowtie -p4 -q -m1 --sam $BIDX $READS/$SAMPLE.umi.fq | samtools view -ShuF4 - | samtools sort -o - $SAMPLE.temp -m 9500000000 > $SAMPLE.has_umi.bam
gzip $READS/$SAMPLE.umi.fq
samtools index $SAMPLE.has_umi.bam

python $HOME/projects/utils/umitools-process-bam.py $SAMPLE.has_umi.bam $SAMPLE.umi_filtered.bam
samtools index $SAMPLE.umi_filtered.bam