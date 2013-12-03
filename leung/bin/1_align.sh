#!/usr/bin/env bash
#BSUB -J align[1-2]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 8
#BSUB -P leung

<<DOC
align rnaseq reads using gsnap
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/leung/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq=$DATA/common/$sample.fastq.gz
results=$RESULTS/$sample
bam=$results/$sample.bam

if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $bam ]]; then
    gsnap -D ~analysiscore/genomes/GMAPDB/hg19_semiTotal -d hg19_semiTotal --gunzip \
        -t 8 -N 1 -s ~analysiscore/genomes/GMAPDB/hg19_semiTotal/hg19_semiTotal.maps/hg19_semiTotal.ensembl.splicesites \
        -v ~analysiscore/genomes/GMAPDB/hg19_semiTotal/hg19_semiTotal.maps/hg19_semiTotal.snp135_strict \
        --read-group-id $sample --read-group-name $sample --read-group-library $sample --read-group-platform Illumina \
        -Q -n 20 -A sam $fastq \
        | samtools view -ShuF4 - \
        | samtools sort -o -m 8G - $sample.temp \
        > $bam
fi
