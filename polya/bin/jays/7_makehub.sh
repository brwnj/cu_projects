#!/usr/bin/env bash
#BSUB -J makehub
#BSUB -e makehub.%J.err
#BSUB -o makehub.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1

<<DOC
pull all necessary files and convert as necessary
DOC

set -o nounset -x
source $HOME/projects/polya/bin/config.sh

cd $HUB/hg19
# copy over BWs
cp $RESULT/*/*.bw .
rm *UMIs_not_removed*
# copy over classified peaks
cp $RESULT/*/*classified.bed.gz .
gunzip -f *.classified.bed.gz
# convert to BB6
for bed in *.classified.bed; do bed2bb.py --type bed6 $CHROM_SIZES $bed; done
# delete classified peaks.bed
rm *.classified.bed


# copy over sites beds
cp $POLYASITES/*sites*.bed.gz .
# convert to BB
gunzip -f *sites*.bed.gz
for bed in *sites*.bed; do bed2bb.py --type bed6 $CHROM_SIZES $bed; done
# these file names must match:
    # MP.sites.c13.bb
    # MP.sites.c1234.bb
    # MP.sites.wholegene.bb
    # PK.sites.c13.bb
    # PK.sites.c1234.bb
    # PK.sites.wholegene.bb
rm *sites*.bed


# this will intentionally fail when there are no new files
# move the dexseq visualization files over
cp $SITESHIFTS/*dexseq.bed .
# convert to BB12
for bed in *dexseq*.bed; do
    bed2bb.py --type bed12 $CHROM_SIZES $bed
done
# delete beds
rm *dexseq*.bed


for f in $FISHERRESULTS/*fisher*txt.gz; do
    fbase=$(basename $f)
    fisherbed=${f/.txt.gz/.bed}
    fisherbb=${f/.txt.gz/.bb}
    if [[ ! -f $fisherbb ]]; then
        # create beds of fisher test results
        python $BIN/visualize_fisher_shifts.py $f $POLYASITES/${fbase:0:2}.sites.c13.bed.gz \
            | bedtools sort -i - > $fisherbed
        # convert bed to bigbed
        bed2bb.py --type bed12 $CHROM_SIZES $fisherbed
        # delete the bed
        rm $fisherbed
    fi
    # cp over the bb
    cp $fisherbb $HUB/hg19/
done


# run generate_trackdb.py > trackDb.txt
python $HOME/projects/polya/bin/generate_trackdb.py . $METADATA > trackDb.txt

# rsync -rvu --delete ~/projects/polya/results/common/hub amc-sandbox:/data/home/brownj/public_html/polya/
