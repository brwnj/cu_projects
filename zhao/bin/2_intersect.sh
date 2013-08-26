#!/usr/bin/env bash
#BSUB -J intersect
#BSUB -e intersect.%J.err
#BSUB -o intersect.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P zhao

<<DOC
process intersection of sups minus wt
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/zhao/bin/config.sh

for f in ACATCG GCCTAA TGGTCA; do
    if [[ ! -f $RESULTS/$f/$f.snpeff.vcf ]]; then
        # record snps from unified genotyper
        java -jar $SNPEFF eff -chr chr -noStats -v -c $SNPEFFCONFIG -snp sacCer3 $RESULTS/$f/$f.ug.vcf > $RESULTS/$f/$f.ug.snpeff.vcf
        # record everything other than snps from haplotype caller
        java -jar $SNPEFF eff -chr chr -noStats -v -c $SNPEFFCONFIG -del -ins sacCer3 $RESULTS/$f/$f.hc.vcf > $RESULTS/$f/$f.hc.snpeff.vcf
        # combine into one vcf
        cp $RESULTS/$f/$f.ug.snpeff.vcf $RESULTS/$f/$f.snpeff.vcf.unsorted
        awk '$1 !~ /#/' $RESULTS/$f/$f.hc.snpeff.vcf >> $RESULTS/$f/$f.snpeff.vcf.unsorted
        bedtools sort -i $RESULTS/$f/$f.snpeff.vcf.unsorted > $RESULTS/$f/$f.snpeff.vcf
        rm $RESULTS/$f/$f.snpeff.vcf.unsorted
    fi
done

bedtools intersect -f .9 -r -a $RESULTS/GCCTAA/GCCTAA.snpeff.vcf -b $RESULTS/TGGTCA/TGGTCA.snpeff.vcf > $RESULTS/non_unique_intersection.vcf
bedtools intersect -f .9 -r -a $RESULTS/TGGTCA/TGGTCA.snpeff.vcf -b $RESULTS/GCCTAA/GCCTAA.snpeff.vcf >> $RESULTS/non_unique_intersection.vcf

uniq $RESULTS/non_unique_intersection.vcf > $RESULTS/intersection.vcf
rm $RESULTS/non_unique_intersection.vcf
bedtools subtract -a $RESULTS/intersection.vcf -b $RESULTS/ACATCG/ACATCG.snpeff.vcf > $RESULTS/intersection_minus_wt.vcf
java -jar ~/opt/snpeff/snpEff.jar eff -chr chr -noStats -v -c ~/opt/snpeff/snpEff.config -o txt sacCer3 $RESULTS/intersection_minus_wt.vcf > $RESULTS/intersection_minus_wt.txt.stupidheader
tail -n +3 $RESULTS/intersection_minus_wt.txt.stupidheader > $RESULTS/intersection_minus_wt.txt
rm $RESULTS/intersection_minus_wt.txt.stupidheader
# python ~/projects/zhao/bin/annotate.py intersection_minus_wt.txt ~/ref/sacCer3/sacCer3.descriptions > intersection_minus_wt.annotated.txt
