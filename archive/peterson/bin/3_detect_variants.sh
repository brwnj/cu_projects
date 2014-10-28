#! /usr/bin/env bash
#BSUB -J detect.variants
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -R "rusage[mem=34] span[hosts=1]"
#BSUB -n 4

set -o nounset -o pipefail -o errexit -x

SAMPLES=(3Poly_1794e3_Mutant)
SAMPLE=${SAMPLES[$LSB_JOBINDEX]}

CORES=4
RAM='-Xmx34g'

# gatk
# 
# realignertargetcreator
# indelrealigner using output from realignertargetcreator
# BaseRecalibrator input using realigned
# #printreads for brent
# unifiedgenotyper

PICARD=$HOME/opt/picard-tools-1.68/picard-tools-1.68
GATK=$HOME/opt/GenomeAnalysisTK-1.6-11-g3b2fab9
REF=/vol1/home/kenjones/genomes/mm9
REFVCF=/vol1/home/gowank/genomes/mm9/mm9_dbSNP.vcf
BIN=/vol1/home/kenjones/bin

BAM=$SAMPLE.bam
NODUPBAM=$SAMPLE.nodup.bam
METRICS=$SAMPLE.metrics
NODUPRG=$SAMPLE.nodup.rg.bam
INTERVALS=$SAMPLE.intervals
REALIGNED=$SAMPLE.nodup.realigned.bam
RECALCSV=$SAMPLE.nodup.realigned.csv
RECALBAM=$SAMPLE.nodup.realigned.recalibrated.bam
NEWRGBAM=$SAMPLE.newrg.bam
VCF=$SAMPLE.combined.snps.indels.raw.vcf
ANNOINPUT=$SAMPLE.combined.snps.indels.avinput


if [ ! -f $NODUPBAM ]; then
    java $RAM -jar $PICARD/MarkDuplicates.jar \
        ASSUME_SORTED=true \
        INPUT=$BAM \
        OUTPUT=$NODUPBAM \
        METRICS_FILE=$METRICS \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=800 \
        CREATE_INDEX=true
fi

# add read groups since gatk does not support without them
if [ ! -f $NODUPRG ]; then
    java $RAM -jar $PICARD/AddOrReplaceReadGroups.jar \
        I=$NODUPBAM \
        O=$NODUPRG \
        SORT_ORDER=coordinate \
        RGID=id \
        RGLB=lb \
        RGPL=illumina \
        RGSM=$SAMPLE \
        RGPU=1 \
        CREATE_INDEX=true
fi

if [ ! -f $INTERVALS ]; then
    java $RAM -jar $GATK/GenomeAnalysisTK.jar \
        --analysis_type RealignerTargetCreator \
        --num_threads $CORES \
        --reference_sequence $REF/mm9.fa \
        --out $INTERVALS \
        --input_file $NODUPRG
fi

if [ ! -f $REALIGNED ]; then
    java $RAM -jar $GATK/GenomeAnalysisTK.jar \
        --analysis_type IndelRealigner \
        --num_threads $CORES \
        --reference_sequence $REF/mm9.fa \
        --input_file $NODUPRG \
        --targetIntervals $INTERVALS \
        --out $REALIGNED
fi

if [ ! -f $RECALCSV ]; then
    java $RAM -jar /vol1/home/kenjones/bin/GenomeAnalysisTK.jar \
    	-l INFO \
    	-R /vol1/home/kenjones/genomes/mm9/mm9.fa \
    	-B:dbsnp,VCF /vol1/home/kenjones/genomes/mm9/mm9_dbSNP.vcf \
    	-I $REALIGNED \
    	--default_platform ILLUMINA \
    	--default_read_group $SAMPLE \
    	-T CountCovariates \
    	-cov ReadGroupCovariate \
    	-cov QualityScoreCovariate \
    	-cov CycleCovariate \
    	-cov DinucCovariate \
    	-recalFile $RECALCSV
        # 
        # java $RAM -jar $GATK/GenomeAnalysisTK.jar \
        #   --analysis_type CountCovariates \
        #   --num_threads $CORES \
        #   --logging_level INFO \
        #   --reference_sequence $REF/mm9.fa \
        #   --knownSites $REFVCF \
        #   --input_file $REALIGNED \
        #   --default_platform ILLUMINA \
        #   --covariate ReadGroupCovariate \
        #   --covariate QualityScoreCovariate \
        #   --covariate CycleCovariate \
        #   --covariate DinucCovariate \
        #   --recal_file $RECALCSV
fi

if [ ! -f $RECALBAM ]; then
    java $RAM -jar $GATK/GenomeAnalysisTK.jar \
    	--analysis_type TableRecalibration \
    	--logging_level INFO \
    	--reference_sequence $REF/mm9.fa \
    	--input_file $REALIGNED \
    	--default_platform ILLUMINA \
    	--out $RECALBAM \
    	--recal_file $RECALCSV
fi

if [ ! -f $NEWRGBAM ]; then
    java $RAM -jar $PICARD/AddOrReplaceReadGroups.jar \
        INPUT=$RECALBAM \
        OUTPUT=$NEWRGBAM \
        SORT_ORDER=coordinate \
        RGID=id \
        RGLB=lb \
        RGPL=illumina \
        RGSM=$SAMPLE \
        RGPU=1 \
        CREATE_INDEX=true
fi

if [ ! -f $VCF ]; then
    # java $RAM -jar $GATK/GenomeAnalysisTK.jar \
    #   --analysis_type UnifiedGenotyper \
    #   --reference_sequence $REF/mm9.fa \
    #   --input_file $NEWRGBAM \
    #   --dbsnp $REFVCF \
    #   --out $VCF \
    #   --genotype_likelihoods_model BOTH \
    #   --standard_min_confidence_threshold_for_calling 50.00 \
    #   --standard_min_confidence_threshold_for_emitting 10.0
    
    java $RAM -jar /vol1/home/kenjones/bin/GenomeAnalysisTK.jar \
    	-R /vol1/home/kenjones/genomes/mm9/mm9.fa \
    	-T UnifiedGenotyper \
    	-I $NEWRGBAM \
    	-B:dbsnp,VCF /vol1/home/kenjones/genomes/mm9/mm9_dbSNP.vcf \
    	-o $VCF \
    	-glm BOTH \
    	-stand_call_conf 50.00 \
    	-stand_emit_conf 10.0
fi

perl $BIN/convert2annovar.pl \
    -format vcf4 \
    -includeinfo \
    -filter pass \
    $VCF > $ANNOINPUT