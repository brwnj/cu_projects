#! /usr/bin/env bash
#BSUB -J novoalign[1-12]
#BSUB -e novoalign.%J.%I.err
#BSUB -o novoalign.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 4

<<DOC
Align using Novoalign, suppressing all reads that align more than once.
DOC

set -o nounset -o pipefail -o errexit -x

# SAMPLES=(idx0 MP51 MP52 MP53 PK61 PK62)
# SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}

samples=(idx0 100-3 86-7 93-1 93-3 m+c m+e2
            nbt29 nbt39 nbt89 ts21 ts28 ts57)
sample=${samples[${LSB_JOBINDEX}]}

fastq=$HOME/projects/polya/data/20130114/$sample.umi.fq.gz
novoidx=$HOME/projects/hits-clip/data/common/novoalign/hg18
umibam=$HOME/projects/polya/results/common/$sample/$sample.UMIs_not_removed.bam
bam=$sample.bam
bin=$HOME/devel/umitools/umitools

novoalign -d $novoidx -f $fastq -o SAM -l 20 -s 5 -r None -c 4 -k \
    2> $sample.alignment.txt \
    | samtools view -ShuF4 - \
    | samtools sort -o - $sample.temp -m 9500000000 \
    > $umibam

samtools index $umibam
python $bin/process_bam.py $umibam $bam
samtools index $bam