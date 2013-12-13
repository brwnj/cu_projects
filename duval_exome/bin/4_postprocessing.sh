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
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/duval_exome/bin/config.sh

intersect_results=$RESULTS/intersect_results
if [[ ! -d $intersect_results ]]; then
    mkdir -p $intersect_results
fi

normal_samples=(172160 353GAAACC 166408)
tumor_samples=(500 353GGCTAC 400)

for (( i = 0; i < ${#normal_samples[@]}; i++ )); do
    normal_sample=${normal_samples[$i]}
    tumor_sample=${tumor_samples[$i]}
    intersect=$intersect_results/$normal_sample.$tumor_sample.vcf
    snpeff=$intersect_results/$normal_sample.$tumor_sample.snpeff.txt
    annotated=$intersect_results/$normal_sample.$tumor_sample.annotated_snps.txt
    if [[ ! -f $intersect ]]; then
        bedtools intersect -v -a $RESULTS/$tumor_sample/$tumor_sample.freebayes.vcf -b $RESULTS/$normal_sample/$normal_sample.freebayes.vcf > $intersect
    fi
    if [[ ! -f $snpeff ]]; then
        java -jar $SNPEFF eff -chr chr -noStats -v -c $SNPEFFCONFIG -o txt -minQ 10 -minC 10 -no-downstream -no-intergenic -no-intron -no-upstream canis_familiaris $intersect > $snpeff
    fi
    if [[ ! -f $annotated ]]; then
        # clip the header
        tmp=$intersect_results/$normal_sample.$tumor_sample.tmp
        tail -n +3 $snpeff > $tmp
        mv $tmp $snpeff
        python $BIN/annotate.py $snpeff $ANNOTATION > $annotated
    fi
done
