#! /usr/bin/env bash
#BSUB -J mirna_seeds
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

MIRBASE=$HOME/projects/hits-clip/data/common/mirbase/mirbase-18/mature.hsa.fa.gz
SEED_START=2
SEED_ENDS="8 9"

for SEED_END in $SEED_ENDS; do
    HSA="hsa.seeds.$SEED_START.$SEED_END.fa"
    cliptools-report-mirna-seeds -s $SEED_START -e $SEED_END $MIRBASE > $HSA
done