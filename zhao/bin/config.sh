#!/usr/bin/env bash
set -o nounset

# wt   ACATCG
# sup6 GCCTAA
# sup8 TGGTCA

# present in sup6 and sup8 and not present in wt

SAMPLES=(ACATCG GCCTAA TGGTCA)
NOVOIDX=$HOME/ref/sacCer3/sacCer3.nix
RESULTS=$HOME/projects/zhao/results/common
DATA=$HOME/projects/zhao/data/common
PICARD=$HOME/opt/picard
GATK=$HOME/opt/gatk/GenomeAnalysisTK.jar
# need both for the reference fasta after sorting by chrom
# samtools faidx $REFERENCE
# java -jar ~/opt/picard/CreateSequenceDictionary.jar R=sacCer3.fa O=sacCer3.dict GENOME_ASSEMBLY=sacCer3 SPECIES=sacCer3
REFERENCE=$HOME/ref/sacCer3/sacCer3.fa

# <<DOC
# 
# for f in ACATCG GCCTAA TGGTCA; do java -jar ~/opt/snpeff/snpEff.jar eff -chr chr -noStats -v -c ~/opt/snpeff/snpEff.config -snp sacCer3 ~/projects/zhao/results/common/$f/$f.ug.vcf > ~/projects/zhao/results/common/$f/$f.ug.snpeff.vcf; done
# for f in ACATCG GCCTAA TGGTCA; do java -jar ~/opt/snpeff/snpEff.jar eff -chr chr -noStats -v -c ~/opt/snpeff/snpEff.config -del -ins sacCer3 ~/projects/zhao/results/common/$f/$f.hc.vcf > ~/projects/zhao/results/common/$f/$f.hc.snpeff.vcf; done
# 
# java -jar ~/opt/snpeff/snpEff.jar eff -chr chr -noStats -v -c ~/opt/snpeff/snpEff.config sacCer3 ~/projects/zhao/results/common/ACATCG/ACATCG.ug.vcf > ~/projects/zhao/results/common/ACATCG/ACATCG.snpeff.vcf
# java -jar ~/opt/snpeff/snpEff.jar eff -chr chr -noStats -v -c ~/opt/snpeff/snpEff.config -o txt sacCer3 ~/projects/zhao/results/common/ACATCG/ACATCG.ug.vcf > ~/projects/zhao/results/common/ACATCG/ACATCG.snpeff.txt
# 
# bedtools intersect -f .9 -a ~/projects/zhao/results/common/GCCTAA/GCCTAA.snpeff.vcf -b ~/projects/zhao/results/common/TGGTCA/TGGTCA.snpeff.vcf > non_unique_intersection.vcf
# bedtools intersect -f .9 -a ~/projects/zhao/results/common/TGGTCA/TGGTCA.snpeff.vcf -b ~/projects/zhao/results/common/GCCTAA/GCCTAA.snpeff.vcf >> non_unique_intersection.vcf
# uniq non_uniq_intersection.vcf > intersection.vcf
# bedtools subtract -a intersection.vcf -b ACATCG/ACATCG.snpeff.vcf > intersection_minus_wt.vcf
# java -jar ~/opt/snpeff/snpEff.jar eff -chr chr -noStats -v -c ~/opt/snpeff/snpEff.config -o txt sacCer3 intersection_minus_wt.vcf > intersection_minus_wt.txt
# DOC