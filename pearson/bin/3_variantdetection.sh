#!/usr/bin/env bash
#BSUB -J variant.detection[1-2]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 1
#BSUB -P pearson

<<DOC
known cause: scf_8253930:65129
Call variants using gatk's haplotype caller
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/pearson/bin/config.sh

name=${NAMES[$(($LSB_JOBINDEX - 1))]}
shortname=$(echo $name | cut -f1 -d_)
java='java -Xmx24g -jar'
outdir=$RESULTS/common/$name
bam=$outdir/$name.bam
rgbam=${bam/.bam/.rdgroup.bam}
nodups=${bam/.bam/.nodups.bam}
duplicatemetrics=${bam/.bam/.dup_metrics.txt}
targetintervals=${bam/.bam/.intervals}
realigned=${bam/.bam/.realign.bam}
vcf=${bam/.bam/.vcf}

# need to add read group info...
# if [ -f $bam ] && [ ! -f $rgbam ]; then
#     $java $PICARD/AddOrReplaceReadGroups.jar \
#         INPUT=$bam \
#         OUTPUT=$rgbam \
#         SORT_ORDER=coordinate \
#         RGID=$(echo $name | cut -f1 -d_) \
#         RGLB=PE \
#         RGPL=illumina \
#         RGPU=$shortname \
#         RGSM=$shortname
# fi
# if [ -f $rgbam ] && [ ! -f $nodups ]; then
#     $java $PICARD/MarkDuplicates.jar \
#         ASSUME_SORTED=true \
#         INPUT=$bam \
#         OUTPUT=$nodups \
#         METRICS_FILE=$duplicatemetrics \
#         CREATE_INDEX=true \
#         MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
#         REMOVE_DUPLICATES=true \
#         MAX_RECORDS_IN_RAM=800000
# fi
# if [ -f $nodups ] && [ ! -f $intervals ]; then
#     $java $GATK --analysis_type RealignerTargetCreator \
#         --reference_sequence $REFERENCE \
#         --input_file $nodups \
#         --out $targetintervals
# fi
# if [ -f $intervals ] && [ ! -f $realigned ]; then
#     $java $GATK --analysis_type IndelRealigner \
#         --reference_sequence $REFERENCE \
#         --input_file $nodups \
#         --targetIntervals $targetintervals \
#         --out $realigned
# fi
if [ -f $realigned ] && [ ! -f $vcf ]; then
    $java $GATK --analysis_type UnifiedGenotyper \
        --input_file $realigned \
        --reference_sequence $REFERENCE \
        --genotype_likelihoods_model BOTH \
        --out $vcf
fi
