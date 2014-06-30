#!/usr/bin/env bash
#BSUB -J align[1-9]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -n 5
#BSUB -q normal

<<DOC
align the reads using rum
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/ribera/bin/ribera.cfg

sample=${SAMPLES[$LSB_JOBINDEX]}
outdir=$RESULTS/common/$sample
if [ ! -d $outdir ]; then
    mkdir -p $outdir
fi

#if [ $CONCATENATE == "TRUE" ]; then
# concatenate reads that fall across lanes
if [ ! -f $DATA/$sample.fastq.gz ] && [ ! -f $DATA/$sample.fastq ]; then
    zcat $DATA/$sample* > $DATA/$sample.fastq
fi

# quality control the reads
if [ ! -f $outdir/${sample}_fastqc/fastqc_report.html ]; then
    $FASTQC --outdir $outdir --threads 4 --format fastq $DATA/$sample.fastq.gz
fi

if [ ! -f $outdir/$sample.bam ]; then
    # rum can't handle gzipped reads
    if [ -f $DATA/$sample.fastq.gz ]; then
        gunzip $DATA/$sample.fastq.gz
    fi
    # align
    rum_runner align -v -i $INDEX -o $outdir --chunks 5 --name $sample $DATA/$sample.fastq

    # cleaning up
    gzip $DATA/$sample.fastq
    cd $outdir
    samtools view -ShuF 4 RUM.sam | samtools sort -o - $sample.temp -m 9500000000 > $sample.bam
    samtools index $sample.bam
    rm -f RUM.sam
    gzip *.fa
    gzip RUM_Unique
    gzip RUM_NU
fi