#! /usr/bin/env bash
#BSUB -J diff
#BSUB -e diff.%J.%I.err
#BSUB -o diff.%J.%I.out
#BSUB -q normal

<<DOC
perform differential expression testing.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/ribera/bin/ribera.cfg

## they need to have common name along with the ensembl id

python ~/projects/utils/cufflinks_anova.py -v -a $RESULTS/common/1A/genes.fpkm_tracking -a $RESULTS/common/2A/genes.fpkm_tracking -a $RESULTS/common/3A/genes.fpkm_tracking \
    -b $RESULTS/common/1B/genes.fpkm_tracking -b $RESULTS/common/2B/genes.fpkm_tracking -b $RESULTS/common/3B/genes.fpkm_tracking > $RESULTS/common/A_vs_B.temp
python ~/projects/utils/qvality.py $RESULTS/common/A_vs_B.temp > $RESULTS/common/A_vs_B.diff
rm $RESULTS/common/A_vs_B.temp

python ~/projects/utils/cufflinks_anova.py -v -a $RESULTS/common/1A/genes.fpkm_tracking -a $RESULTS/common/2A/genes.fpkm_tracking -a $RESULTS/common/3A/genes.fpkm_tracking \
    -b $RESULTS/common/1C/genes.fpkm_tracking -b $RESULTS/common/2C/genes.fpkm_tracking -b $RESULTS/common/3C/genes.fpkm_tracking > $RESULTS/common/A_vs_C.temp
python ~/projects/utils/qvality.py $RESULTS/common/A_vs_C.temp > $RESULTS/common/A_vs_C.diff
rm $RESULTS/common/A_vs_C.temp

python ~/projects/utils/cufflinks_anova.py -v -a $RESULTS/common/1B/genes.fpkm_tracking -a $RESULTS/common/2B/genes.fpkm_tracking -a $RESULTS/common/3B/genes.fpkm_tracking \
    -b $RESULTS/common/1C/genes.fpkm_tracking -b $RESULTS/common/2C/genes.fpkm_tracking -b $RESULTS/common/3C/genes.fpkm_tracking > $RESULTS/common/B_vs_C.temp
python ~/projects/utils/qvality.py $RESULTS/common/B_vs_C.temp > $RESULTS/common/B_vs_C.diff
rm $RESULTS/common/B_vs_C.temp