#!/usr/bin/env bash
#BSUB -J post[1-16]
#BSUB -e post.%J.%I.err
#BSUB -o post.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>1] rusage[mem=1] span[hosts=1]"
#BSUB -n 1
#BSUB -P duval_exome

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/duval_exome/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
snpeff=$RESULTS/$sample/$sample.snpeff.txt.gz

for c in 100 50 20 10; do
    result=$RESULTS/$sample/${sample}_nonsyn_snp_counts_${c}.txt
    python ~/projects/duval_exome/bin/gene_snp_counts.py $snpeff ~/ref/canFam3/ensemblToGeneName.txt.gz --coverage $c > $result
done
