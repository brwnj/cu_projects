#!/usr/bin/env bash
#BSUB -J duval_variants[1-16]
#BSUB -e duval_variants.%J.%I.err
#BSUB -o duval_variants.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 1
#BSUB -P duval_exome

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/duval_exome/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

java='java -Xmx24g -jar'
bam=$RESULTS/$sample/$sample.bam
rgbam=${bam/.bam/.rg.bam}
nodups=${bam/.bam/.nodups.bam}
reordered=${bam/.bam/.reordered.bam}
duplicatemetrics=${bam/.bam/.dup_metrics.txt}
targetintervals=${bam/.bam/.intervals}
realigned=${bam/.bam/.realign.bam}
vcf=${bam/.bam/.vcf}
freebayesvcf=${bam/.bam/.freebayes.vcf}

# if [[ -f $vcf ]]; then
#     echo "processing complete for $sample"
#     exit 0
# fi

# picard
READ_GROUP_OPTS="INPUT=$bam OUTPUT=$rgbam SORT_ORDER=coordinate RGID=$sample RGLB=PE RGPL=illumina RGPU=$sample RGSM=$sample"
MARK_DUPLICATES_OPTS="ASSUME_SORTED=true INPUT=$rgbam OUTPUT=$nodups METRICS_FILE=$duplicatemetrics CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=true MAX_RECORDS_IN_RAM=800000"
REORDER_SAM_OPTS="INPUT=$nodups OUTPUT=$reordered REFERENCE=$REFERENCE"

# gatk
REALIGNER_TARGET_CREATOR="--analysis_type RealignerTargetCreator --reference_sequence $REFERENCE --input_file $reordered --out $targetintervals"
INDEL_REALIGNER="--analysis_type IndelRealigner --reference_sequence $REFERENCE --input_file $reordered --targetIntervals $targetintervals --out $realigned"
UNIFIED_GENOTYPER="--analysis_type UnifiedGenotyper --input_file $realigned --reference_sequence $REFERENCE --genotype_likelihoods_model BOTH --out $vcf"

if [ -f $bam ] && [ ! -f $rgbam ]; then
    $java $PICARD/AddOrReplaceReadGroups.jar $READ_GROUP_OPTS
fi
if [ -f $rgbam ] && [ ! -f $nodups ]; then
    $java $PICARD/MarkDuplicates.jar $MARK_DUPLICATES_OPTS
fi
# reorder the bam to match the reference contig order
if [ -f $nodups ] && [ ! -f $reordered ]; then
    $java $PICARD/ReorderSam.jar $REORDER_SAM_OPTS
    # must index the new bam
    samtools index $reordered
fi
if [ -f $reordered ] && [ ! -f $targetintervals ]; then
    $java $GATK $REALIGNER_TARGET_CREATOR
fi
if [ -f $targetintervals ] && [ ! -f $realigned ]; then
    $java $GATK $INDEL_REALIGNER
fi
if [ -f $realigned ] && [ ! -f $vcf ]; then
    $java $GATK $UNIFIED_GENOTYPER
fi
if [ -f $realigned ] && [ ! -f $freebayesvcf ]; then
    freebayes -b $realigned -v $freebayesvcf -f $REFERENCE -0
fi
