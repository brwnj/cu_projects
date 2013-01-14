#! /usr/bin/env bash
#BSUB -J bowtie[1-5]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

<<DOC
Align using Novoalign, suppressing all reads that align more than once.
DOC

set -o nounset -o pipefail -o errexit -x

# SAMPLES=(idx0 MP51 MP52 MP53 PK61 PK62)
# SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}

samples=()
sample=${samples[${LSB_JOBINDEX}]}

fastq=$HOME/projects/polya/data/??????????

novoidx=$HOME/projects/hits-clip/data/common/novoalign/hg18
umibam
bam

 | samtools view -ShuF4 - | samtools sort -o - $SAMPLE.temp -m 9500000000 > $SAMPLE.has_umi.bam

samtools index $umibam
python $HOME/projects/utils/umitools-process-bam.py $umibam $bam
samtools index $bam