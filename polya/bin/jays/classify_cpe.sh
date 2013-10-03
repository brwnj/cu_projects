#!/usr/bin/env bash
#BSUB -J classify_cpe
#BSUB -e classify_cpe.%J.err
#BSUB -o classify_cpe.%J.out
#BSUB -q normal

<<DOC
Determie CPE elements for the genome
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/polya/bin/config.sh

# using polya sites as reference positions
BED=$RESULT/polya_sites/PC.sites.c13.bed.gz
OUTBED=$RESULT/polya_sites/cpe.sites.c13.bed
python $BIN/classify_cpe.py $BED $FASTA --verbose \
    > $OUTBED

# select unique lines, overlapping genes cause redundant annotations
tmpbed=$(mktemp)
sort $OUTBED | uniq > $tmpbed
mv $tmpbed $OUTBED

# sort for big bed creation
bedSort $OUTBED $OUTBED

bed2bb.py --type "bed9" $CHROM_SIZES $OUTBED
gzip -f $OUTBED
