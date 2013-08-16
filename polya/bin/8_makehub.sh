#!/usr/bin/env bash
#BSUB -J makehub
#BSUB -e makehub.%J.err
#BSUB -o makehub.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

<<DOC
pull all necessary files and convert as necessary
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

cd $HUB/hg19
# copy over BWs
cp $RESULT/*/*.bw .
rm *UMIs_not_removed*
# copy over classified peaks
cp $RESULT/*/*classified.bed.gz .
gunzip *.classified.bed.gz
# convert to BB6
for bed in *.classified.bed; do bed2bb.py --type bed6 $CHROM_SIZES $bed; done
# delete classified peaks.bed
rm *.classified.bed
# copy over sites beds
cp $POLYASITES/*sites*.bed.gz .
# convert to BB
gunzip *sites*.bed.gz
for bed in *sites*.bed; do bed2bb.py --type bed6 $CHROM_SIZES $bed; done
# these file names must match:
    # MP.sites.c13.bb
    # MP.sites.c1234.bb
    # MP.sites.wholegene.bb
    # PK.sites.c13.bb
    # PK.sites.c1234.bb
    # PK.sites.wholegene.bb
rm *sites*.bed
# move the dexseq visualization files over
mv $SITESHIFTS/*.bed .
# convert to BB12
for bed in *_to_*.bed; do bed2bb.py --type bed12 $CHROM_SIZES $bed; done
# delete beds
rm *_to_*.bed
# run generate_trackdb.py > trackDb.txt
python $HOME/projects/polya/bin/generate_trackdb.py . $METADATA > trackDb.txt
