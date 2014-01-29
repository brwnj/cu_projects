#!/usr/bin/env bash
#BSUB -J mirna_abundance
#BSUB -e mirna_abundance.%J.err
#BSUB -o mirna_abundance.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P nicoli

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

python ~/projects/nicoli/bin/mirna_abundance.py ~/ref/mirbase/20/mature.hsa.bed.gz $RESULTS/L1/L1.rmd_pos.seed8_mapping.txt.gz $RESULTS/L1/L1.rmd_neg.seed8_mapping.txt.gz $GENOMEDATA L1_pos L1_neg | sort -k1,1 -k2,2n > $RESULTS/L1/L1_mirna_abundance.bed
python ~/projects/nicoli/bin/mirna_abundance.py ~/ref/mirbase/20/mature.hsa.bed.gz $RESULTS/L1/L1.rmd_pos.seed8_mapping.txt.gz $RESULTS/L1/L1.rmd_neg.seed8_mapping.txt.gz $GENOMEDATA L1.rmd_pos L1.rmd_neg | sort -k1,1 -k2,2n > $RESULTS/L1/L1_rmd_mirna_abundance.bed

python ~/projects/nicoli/bin/mirna_abundance.py ~/ref/mirbase/20/mature.hsa.bed.gz $RESULTS/L2/L2.rmd_pos.seed8_mapping.txt.gz $RESULTS/L2/L2.rmd_neg.seed8_mapping.txt.gz $GENOMEDATA L2_pos L2_neg | sort -k1,1 -k2,2n > $RESULTS/L2/L2_mirna_abundance.bed
python ~/projects/nicoli/bin/mirna_abundance.py ~/ref/mirbase/20/mature.hsa.bed.gz $RESULTS/L2/L2.rmd_pos.seed8_mapping.txt.gz $RESULTS/L2/L2.rmd_neg.seed8_mapping.txt.gz $GENOMEDATA L2.rmd_pos L2.rmd_neg | sort -k1,1 -k2,2n > $RESULTS/L2/L2_rmd_mirna_abundance.bed
