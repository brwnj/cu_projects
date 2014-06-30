#!/usr/bin/env bash
#BSUB -J align[1-24]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 8
#BSUB -P walter_miRNA

<<DOC
align using novoalign
    trim truseq adapter
    minimum output quality as 30
    trim back to 18 bp attempting to align
    report all alignments
    homopolymer filter as 60
    miRNA paired alignment scoring (ZH flag)
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/walter/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq=$DATA/$sample.fastq.gz

if [[ ! -d $RESULTS/$sample ]]; then
    mkdir -p $RESULTS/$sample
fi
# serialized alignment using novoalign
for strain in H37Rv HN878; do
    idx=$HOME/ref/tuberculosis/$strain.novoidx
    # align using hairpin option; simply adds an additional score.
    bam=$RESULTS/$sample/$sample.$strain.bam
    align_metrics=$sample.$strain.alignment.txt
    if [[ ! -f $bam ]]; then
        # keep all alignments expecting miRNAs to align to multiple places
        # will be able to filter some noise using ZH (hairpin score) in sam
        novoalign -d $idx -f $fastq -a CTGGAATTCTCGGGTGCCAAGG -o SAM -t 30 -s 3 -l 18 -r A -c 8 -h 60 -m -k \
            2> $align_metrics \
            | samtools view -ShuF4 - \
            | samtools sort -o -m 8G - $sample.$strain.temp \
            > $bam
    fi
done
