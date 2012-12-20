#!/usr/bin/env bash
#BSUB -J variant.detection[1-2]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 1

<<DOC
Call variants using gatk's haplotype caller

need the fasta reference sorted by contig name
bioawk -c fastx '{print}' in.fa | sort -k1,1n | awk '{print ">"$1;print $2}'
then you need to make it a valid fasta
java -jar picard-tools-1.79/NormalizeFasta.jar I=tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.fasta  O=tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.normalized.fasta
index it
samtools faidx tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.normalized.fasta
create a dict of the fasta
java -jar picard-tools-1.79/CreateSequenceDictionary.jar R=tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.normalized.fasta GENOME_ASSEMBLY=tta_mic O=tetrahymena_thermophila_sb210__mic__2_supercontigs.sorted.normalized.fasta.dict

once complete
cd /vol1/home/brownj/projects/pearson/results/common
bedtools subtract -a DisA1_TGACCA_L001/DisA1_TGACCA_L001.mic.vcf -b B1868_ATCACG_L001/B1868_ATCACG_L001.mic.vcf > DisA1.uniq.mic.vcf

snpEff reference
add entry to snpEff.config
mv gtf to ~/ref/snpeff/<species>
mv the fasta annotation to genomes in ~/ref/snpeff/genomes
java -jar ~/opt/snpeff/snpEff.jar build -gtf22 -v -c ~/opt/snpeff/snpEff.config tta_mic

run snpEff
java -jar ~/opt/snpeff/snpEff.jar -i vcf -o txt -chr supercontig_ -minC 10 -minQ 30 -no-intergenic -ud 200 -hom -c ~/opt/snpeff/snpEff.config tta_mic DisA1.uniq.mic.vcf | awk '$13!="INTRON"' > DisA1_tta_mic.txt
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/pearson/bin/base.cfg

name=${NAMES[$LSB_JOBINDEX]}
java='java -Xmx24g -jar'
outdir=$RESULTS/common/$name
bam=$outdir/$name.mic.bam
# mark duplicates
nodups=${bam/.bam/.nodups.bam}
duplicatemetrics=${bam/.bam/.metrics.txt}
# realignertargetcreator
intervals=${bam/.bam/.intervals}
# indelrealigner
realigned=${bam/.bam/.realign.bam}
# baserecalibrator
grp=${bam/.bam/.grp}
# haplotypecaller
vcf=${bam/.bam/.vcf}
# annovar input conversion
annovar=${bam/.bam/.annovar}

if [ -f $bam ] && [ ! -f $nodups ]; then
    $java $PICARD/MarkDuplicates.jar ASSUME_SORTED=true INPUT=$bam OUTPUT=$nodups \
        METRICS_FILE=$duplicatemetrics CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 REMOVE_DUPLICATES=true MAX_RECORDS_IN_RAM=800000
fi
if [ -f $nodups ] && [ ! -f $intervals ]; then
    $java $GATK --analysis_type RealignerTargetCreator \
        --reference_sequence $REFERENCE --input_file $nodups --out $intervals
fi
if [ -f $intervals ] && [ ! -f $realigned ]; then
    $java $GATK --analysis_type IndelRealigner --reference_sequence $REFERENCE \
        --input_file $nodups --targetIntervals $intervals --out $realigned
fi
# i don't have a vcf for this species
# if [ -f $realigned ] && [ ! -f $grp ]; then
#     $java $GATK --analysis_type BaseRecalibrator --input_file $realigned \
#         --reference_sequence $REFERENCE --out $grp
# fi
if [ ! -f $vcf ]; then # [ -f $grp ] && --BQSR $grp \
    $java $GATK --analysis_type UnifiedGenotyper --intervals $intervals --input_file $realigned \
        --reference_sequence $REFERENCE --validation_strictness STRICT \
        --out $vcf -glm BOTH -stand_call_conf 50.00 -stand_emit_conf 10.0
fi