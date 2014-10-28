#!/usr/bin/env bash
#BSUB -J align[1-6]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P pillai

<<DOC
align rna-seq data using gsnap
DOC

set -o nounset -o pipefail -o errexit -x

samples=(MP64 MP65 MP66 MP67 MP68 MP69)
sample=${samples[$(($LSB_JOBINDEX - 1))]}

reads=$HOME/projects/mp_rnaseq/data
fastq=$reads/$sample.fastq.gz
results=$HOME/projects/mp_rnaseq/results/common/$sample
if [[ ! -d $results ]]; then
    mkdir -p $results
fi
gmapdb=/vol2/home/analysiscore/genomes/GMAPDB
splicesites=hg18.refseq.splicesites
bam=$results/$sample.bam
if [[ ! -f $bam ]]; then
    gsnap -D $gmapdb -d hg18 --gunzip -s $splicesites -n1 -B5 -Q --nofails \
        -t8 -A sam --read-group-id=$sample --read-group-name=$sample $fastq \
        | samtools view -ShuF 4 - \
        | samtools sort -o - $sample.temp -m 8G > $bam
    samtools index $bam
fi
