#! /usr/bin/env bash
#BSUB -J bowtie2[1-20]
#BSUB -e bowtie2.%J.%I.err
#BSUB -o bowtie2.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 4

<<DOC
align the reads using bowtie 2.
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(workaround
10_AAGCTA_L006_R1_001.fastq.gz
11_GTAGCC_L006_R1_001.fastq.gz
12_TACAAG_L006_R1_001.fastq.gz
13_TTGACT_L006_R1_001.fastq.gz
14_GGAACT_L006_R1_001.fastq.gz
15_TGACAT_L006_R1_001.fastq.gz
16_GGACGG_L006_R1_001.fastq.gz
17_GCGGAC_L006_R1_001.fastq.gz
18_TTTCAC_L006_R1_001.fastq.gz
19_GGCCAC_L006_R1_001.fastq.gz
1_CGTGAT_L006_R1_001.fastq.gz
20_CGAAAC_L006_R1_001.fastq.gz
21_TGGTCA_L006_R1_001.fastq.gz
22_TCAAGT_L006_R1_001.fastq.gz
23_ACATCG_L006_R1_001.fastq.gz
24_ATTGGC_L006_R1_001.fastq.gz
3_GCCTAA_L006_R1_001.fastq.gz
5_CACTGT_L006_R1_001.fastq.gz
7_GATCTG_L006_R1_001.fastq.gz
9_CTGATC_L006_R1_001.fastq.gz
)

SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}
BIDX=$HOME/projects/ref/hg19/bowtie2/hg19
NAME=$(basename $SAMPLE .fastq.gz)
FASTQ=$(basename $SAMPLE .gz)
PROJECT=$HOME/projects/thaikoottathil

if [ ! -f $NAME.bam ]; then
    gunzip $PROJECT/data/$SAMPLE
    
    # gatk requires read group to be set
    bowtie2 -p4 -x $BIDX -q --very-sensitive --rg-id $NAME \
        $PROJECT/data/$FASTQ \
        | samtools view -bSh - > $NAME.unsorted.bam

    samtools sort $NAME.unsorted.bam $NAME

    samtools index $NAME.bam

    gzip $PROJECT/data/$FASTQ
fi