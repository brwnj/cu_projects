#! /usr/bin/env bash
#BSUB -J bowtie2.align[1-7]
#BSUB -e bowtie2.%J.%I.err
#BSUB -o bowtie2.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>21] rusage[mem=20] span[hosts=1]"
#BSUB -n 4

<<DOC
align fastqs using bowtie then sort and index the output. it also generates
bigwigs for use in UCSC.
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(workaround
1_TCAG_L007_R1_001.fastq
2_AGCA_L007_R1_001.fastq
3_CTGT_L007_R1_001.fastq
4_GATC_L007_R1_001.fastq
5_TGAG_L007_R1_001.fastq
6_ACTC_L007_R1_001.fastq
7_CAGT_L007_R1_001.fastq
)

SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}
BIDX=$HOME/projects/ref/sacCer3/sacCer3
NAME=$(basename $SAMPLE .fastq)
PROJECT=$HOME/projects/xxxxxxxx

if [ ! -f $NAME.bam ]; then
    gunzip $PROJECT/data/20120615/$SAMPLE.gz
    
    # gatk requires read group to be set
    bowtie2 -p4 -x $BIDX -q --very-sensitive --rg-id $NAME \
        $PROJECT/data/20120615/$SAMPLE \
        | samtools view -bSh - > $NAME.unsorted.bam

    samtools sort $NAME.unsorted.bam $NAME
    rm $NAME.unsorted.bam
    samtools index $NAME.bam
    gzip $PROJECT/data/20120615/$SAMPLE
fi