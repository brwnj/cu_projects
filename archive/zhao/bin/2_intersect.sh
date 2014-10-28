#!/usr/bin/env bash
#BSUB -J intersect
#BSUB -e intersect.%J.err
#BSUB -o intersect.%J.out
#BSUB -q short
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P zhao

<<DOC
process intersection of sups minus wt
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/zhao/bin/config.sh

wt=$RESULTS/ACATCG/ACATCG.freebayes.vcf
sup6=$RESULTS/GCCTAA/GCCTAA.freebayes.vcf
sup8=$RESULTS/TGGTCA/TGGTCA.freebayes.vcf

sup6sub=${sup6/.vcf/.remaining.vcf}
sup8sub=${sup8/.vcf/.remaining.vcf}

sup6badhead=${sup6/.vcf/.badheader.txt}
sup8badhead=${sup8/.vcf/.badheader.txt}

sup6preanno=${sup6/.vcf/.snpeff.txt}
sup8preanno=${sup8/.vcf/.snpeff.txt}

sup6anno=${sup6/.vcf/.annotated.txt}
sup8anno=${sup8/.vcf/.annotated.txt}

if [[ ! -f $sup6sub ]] && [[ ! -f $sup8sub ]]; then
    # subtract wt from sup samples individually
    bedtools subtract -A -a $sup6 -b $wt > $sup6sub
    bedtools subtract -A -a $sup8 -b $wt > $sup8sub
fi

if [[ ! -f $sup6preanno ]] && [[ ! -f $sup8preanno ]]; then
    # annotate remaining mutations with snpeff
    java -jar $SNPEFF eff -chr chr -noStats -v -c $SNPEFFCONFIG -o txt sacCer3 $sup6sub > $sup6badhead
    java -jar $SNPEFF eff -chr chr -noStats -v -c $SNPEFFCONFIG -o txt sacCer3 $sup8sub > $sup8badhead
    tail -n +3 $sup6badhead > $sup6preanno
    tail -n +3 $sup8badhead > $sup8preanno
    rm -f $sup6badhead $sup8badhead
fi

# add detailed descriptions to both
python $BIN/annotate.py $sup6preanno $ANNOTATION > $sup6anno
python $BIN/annotate.py $sup8preanno $ANNOTATION > $sup8anno

# combine with sample id column
# awk '{split(FILENAME, fn, "/"); print fn[1], $0}' $sup6anno $sup8anno > $RESULTS/annotated_hc_variants.txt
# need to fix second header being added
# need to sort
