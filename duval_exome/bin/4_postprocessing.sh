#!/usr/bin/env bash
#BSUB -J postprocessing
#BSUB -e postprocessing.%J.err
#BSUB -o postprocessing.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P duval_exome

<<DOC
overlap the cancer samples and filter down to most putative mutations

Tumor   Normal
500     172160
353     353M2
400     166408

Sample 353 is tagged GGCTAC
Sample 353M2 is tagged GAAACC
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/duval_exome/bin/config.sh

intersect_results=$RESULTS/intersect_results
if [[ ! -d $intersect_results ]]; then
    mkdir -p $intersect_results
fi

normal_samples=(172160 353GAAACC 166408)
# the first three are paired with normals
tumor_samples=(500 353GGCTAC 400 1025 113 22 36 522 730 735 868 Bililey KTCC)

# run normal/cancer pairs with filtering
for (( i = 0; i < ${#normal_samples[@]}; i++ )); do
    normal=${normal_samples[$i]}
    tumor=${tumor_samples[$i]}
    vcf=$intersect_results/filtering_with_pair_only/$tumor.vcf
    snpeff=$intersect_results/filtering_with_pair_only/$tumor.snpeff.txt.gz
    annotated=$intersect_results/filtering_with_pair_only/$tumor.annotated_snps.txt.gz
    if [[ ! -f $vcf ]]; then
        bedtools intersect -v -a $RESULTS/$tumor/$tumor.freebayes.vcf -b $RESULTS/$normal/$normal.freebayes.vcf > $vcf
    fi
    if [[ ! -f $snpeff ]]; then
        java -jar $SNPEFF eff -chr chr -noStats -v -c $SNPEFFCONFIG \
            -fi $TARGETS -o txt -minQ 10 -minC 10 -no-downstream -no-intergenic \
            -no-intron -no-upstream canis_familiaris $vcf | gzip -c > $snpeff
    fi
    if [[ ! -f $annotated ]]; then
        python $BIN/annotate.py $snpeff $ANNOTATION | gzip -c > $annotated
    fi
done

# create a vcf of normal intersection
filtered_vcfs=""
for (( i = 0; i < ${#normal_samples[@]}; i++ )); do
    normal=${normal_samples[$i]}
    vcf=$RESULTS/$normal/$normal.freebayes.vcf
    filtered=$RESULTS/$normal/$normal.filtered.freebayes.vcf
    filtered_vcfs=$filtered" "$filtered_vcfs
    if [[ ! -f $filtered ]]; then
        cat $vcf | java -jar $SNPSIFT filter "((QUAL >= 500) & (DP >= 100))" | java -jar $SNPSIFT intervals $TARGETS > $filtered
    fi
done

# combine the normal VCFs into a single
if [[ ! -f $NORMALVCF ]]; then
    cat $filtered_vcfs >> tmp.vcf
    bedtools sort -i tmp.vcf > $NORMALVCF
    rm tmp.vcf
fi

# run all samples against normal database
for (( i = 0; i < ${#tumor_samples[@]}; i++ )); do
    tumor=${tumor_samples[$i]}
    vcf=$intersect_results/$tumor.vcf
    snpeff=$intersect_results/$tumor.snpeff.txt.gz
    annotated=$intersect_results/$tumor.annotated_snps.txt.gz
    if [[ ! -f $vcf ]]; then
        bedtools intersect -v -a $RESULTS/$tumor/$tumor.freebayes.vcf -b $NORMALVCF > $vcf
    fi
    if [[ ! -f $snpeff ]]; then
        java -jar $SNPEFF eff -noStats -v -c $SNPEFFCONFIG \
            -fi $TARGETS -o txt -minQ 100 -minC 10 -no-downstream -no-intergenic \
            -no-intron -no-upstream canis_familiaris $vcf | gzip -c > $snpeff
    fi
    if [[ ! -f $annotated ]]; then
        python $BIN/annotate.py $snpeff $ANNOTATION | gzip -c > $annotated
    fi
done
