#!/usr/bin/env bash
#BSUB -J makehub
#BSUB -e makehub.%J.err
#BSUB -o makehub.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

<<DOC
pull all necessary files and convert everything as necessary
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

cd $HUB
# copy over BWs
cp $RESULT/*/*.bw .
rm *UMIs_not_removed*
# copy over classified peaks
cp $RESULT/*/*classified.bed.gz
# convert to BB6

# delete classified peaks.bed

# copy over sites beds
# convert to BB
# ensure names are 
MP.sites.c13.bb
MP.sites.c1234.bb
MP.sites.wholegene.bb
PK.sites.c13.bb
PK.sites.c1234.bb
PK.sites.wholegene.bb

# copy over dexseq_shifts
# convert to BB12
# delete beds

# run generate_trackdb.py > trackDb.txt
