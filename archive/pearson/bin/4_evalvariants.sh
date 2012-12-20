#!/usr/bin/env bash
#BSUB -J gatk.evalvariants
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 1

<<DOC
compare vcfs from unified genotyper
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/pearson/bin/base.cfg

java='java -Xmx24g -jar'
outdir=$RESULTS/common/var_comp.txt
evalvcf=$RESULTS/common/B1868_ATCACG_L001/B1868_ATCACG_L001.vcf
compvcf=$RESULTS/common/DisA1_TGACCA_L001/DisA1_TGACCA_L001.vcf

$java $GATK --analysis_type VariantEval --eval $evalvcf --comp $compvcf --reference_sequence $REFERENCE --out $outdir