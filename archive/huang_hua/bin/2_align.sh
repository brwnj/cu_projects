#!/usr/bin/env bash
#BSUB -J align[1-6]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -n 5
#BSUB -q normal

<<DOC
align the reads using rum
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/huang_hua/bin/huang_hua.cfg

sample=${SAMPLES[$LSB_JOBINDEX]}
outdir=$RESULTS/common/$sample

if [ ! -f ${outdir}/RUM.sam ]; then
    gunzip $DATA/$sample.fastq.gz
    rum_runner align -v -i $INDEX -o $outdir --chunks 5 --name $sample $DATA/$sample.fastq
    gzip $DATA/$sample.fastq
fi

# cleaning up
cd $outdir
samtools view -ShuF 4 RUM.sam | samtools sort -o - $sample.temp -m 9500000000 > $sample.bam
samtools index $sample.bam
rm -f RUM.sam
gzip *.fa
gzip RUM_Unique
gzip RUM_NU
gzip *junction*