#!/usr/bin/env bash
#BSUB -J detect_vars[1-3]
#BSUB -e detect_vars.%J.%I.err
#BSUB -o detect_vars.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 1
#BSUB -P zhao

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/zhao/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
java='java -Xmx24g -jar'
outdir=$RESULTS/common/$sample
bam=$outdir/$sample.bam
nodups=${bam/.bam/.nodups.bam}
duplicatemetrics=${bam/.bam/.dup_metrics.txt}
targetintervals=${bam/.bam/.intervals}
realigned=${bam/.bam/.realign.bam}
vcf=${bam/.bam/.vcf}

if [ ! -f $nodups ]; then
    $java $PICARD/MarkDuplicates.jar \
        ASSUME_SORTED=true \
        INPUT=$bam \
        OUTPUT=$nodups \
        METRICS_FILE=$duplicatemetrics \
        CREATE_INDEX=true \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        REMOVE_DUPLICATES=true \
        MAX_RECORDS_IN_RAM=800000
fi
if [ -f $nodups ] && [ ! -f $intervals ]; then
    $java $GATK --analysis_type RealignerTargetCreator \
        --reference_sequence $REFERENCE \
        --input_file $nodups \
        --out $targetintervals
fi
if [ -f $intervals ] && [ ! -f $realigned ]; then
    $java $GATK --analysis_type IndelRealigner \
        --reference_sequence $REFERENCE \
        --input_file $nodups \
        --targetIntervals $targetintervals \
        --out $realigned
fi
if [ -f $realigned ] && [ ! -f $vcf ]; then
    $java $GATK --analysis_type UnifiedGenotyper \
        --input_file $realigned \
        --reference_sequence $REFERENCE \
        --genotype_likelihoods_model BOTH \
        --out $vcf
fi
