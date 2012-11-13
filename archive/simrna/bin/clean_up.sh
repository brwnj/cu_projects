#! /usr/bin/env bash
#BSUB -J toomuchrum.out[1-8]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

<<DOC
clean up the rum output
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(idx0
Human_RNAseq_gene_rep_1_sample_1_shred.fastq
Human_RNAseq_gene_rep_1_sample_2_shred.fastq
Human_RNAseq_gene_rep_1_sample_3_shred.fastq
Human_RNAseq_gene_rep_1_sample_4_shred.fastq
Human_RNAseq_gene_rep_1_sample_5_shred.fastq
Human_RNAseq_gene_rep_1_sample_6_shred.fastq
Human_RNAseq_gene_rep_1_sample_7_shred.fastq
Human_RNAseq_gene_rep_1_sample_8_shred.fastq
)
SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}
NAME=$(basename $SAMPLE .fastq)
OUTPUT=$HOME/projects/simrna/results/common/$NAME

# cleaning up
cd $OUTPUT
samtools view -Su RUM.sam | samtools sort -o - deleteme > RUM.bam
samtools index RUM.bam
rm -f RUM.sam
gzip *.fa
gzip RUM_Unique
gzip RUM_NU