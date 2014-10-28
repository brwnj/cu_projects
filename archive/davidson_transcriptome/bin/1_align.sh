#!/usr/bin/env bash
#BSUB -J align[1-10]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P davidson_transcriptome

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/davidson_transcriptome/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq=$DATA/$sample.fastq.gz
results=$RESULTS/$sample
bam=$results/$sample.bam

if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $bam ]]; then
    # tophat2 -o $results -p 8 -G $GENES --transcriptome-index $TRANSCRIPTOME $BTBASE $fastq
    # mv $results/accepted_hits.bam $bam
    gsnap -D ~analysiscore/genomes/GMAPDB/hg19_total -d hg19_total --gunzip \
        -t 8 -N 1 -s ~analysiscore/genomes/GMAPDB/hg19_total/hg19_total.maps/hg19_total.ensembl.splicesites \
        -Q -n 20 -A sam $fastq \
        | samtools view -ShuF4 - \
        | samtools sort -o -m 8G - $sample.temp \
        > $bam
fi
