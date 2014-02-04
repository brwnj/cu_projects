#!/usr/bin/env bash
#BSUB -J variants[1-24]
#BSUB -e variants.%J.%I.err
#BSUB -o variants.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P leung

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/leung/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

java='java -Xmx16g -jar'
base_name=$RESULTS/$sample/$sample
bam=$base_name.bam
rgbam=$base_name.rg.bam
nodups=$base_name.nodups.bam
reordered=$base_name.reordered.bam
duplicatemetrics=$base_name.dup_metrics.txt
vcf=$base_name.vcf

if [[ -f $vcf ]]; then
    echo "processing complete for $sample"
    exit 0
fi

# picard
MARK_DUPLICATES_OPTS="ASSUME_SORTED=true INPUT=$rgbam OUTPUT=$nodups METRICS_FILE=$duplicatemetrics CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=true MAX_RECORDS_IN_RAM=800000"
REORDER_SAM_OPTS="INPUT=$nodups OUTPUT=$reordered REFERENCE=$REFERENCE"

if [ -f $rgbam ] && [ ! -f $nodups ]; then
    $java $PICARD/MarkDuplicates.jar $MARK_DUPLICATES_OPTS
fi

# reorder the bam to match the reference contig order
if [ -f $nodups ] && [ ! -f $reordered ]; then
    $java $PICARD/ReorderSam.jar $REORDER_SAM_OPTS
    # must index the new bam
    samtools index $reordered
fi

if [ -f $reordered ] && [ ! -f $vcf ]; then
    freebayes -b $realigned -v $vcf -f $REFERENCE -0
fi
