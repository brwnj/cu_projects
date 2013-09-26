#!/usr/bin/env bash
#BSUB -J run_pooled
#BSUB -e run_pooled.%J.err
#BSUB -o run_pooled.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

<<DOC
i was just running this script within a results dir and outputting files locally
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

# pool specified counts
python $HOME/projects/polya/bin/get_pool_counts.py $METADATA $HOME/projects/polya/results/common/*/*.counts.txt.gz

# run fisher tests across desired comparisons
bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.neg.txt.gz TS-ERP.neg.txt.gz | gzip -c > NBT_to_TS-ERP.neg.fisher.txt.gz" &
bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.pos.txt.gz TS-ERP.pos.txt.gz | gzip -c > NBT_to_TS-ERP.pos.fisher.txt.gz" &
bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.neg.txt.gz TS.neg.txt.gz | gzip -c > NBT_to_TS.neg.fisher.txt.gz" &
bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.pos.txt.gz TS.pos.txt.gz | gzip -c > NBT_to_TS.pos.fisher.txt.gz" &
bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.neg.txt.gz TS-ERN.neg.txt.gz | gzip -c > NBT_to_TS-ERN.neg.fisher.txt.gz" &
bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py NBT.pos.txt.gz TS-ERN.pos.txt.gz | gzip -c > NBT_to_TS-ERN.pos.fisher.txt.gz" &
bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py TS-ERP.neg.txt.gz TS-ERN.neg.txt.gz | gzip -c > TS-ERP_to_TS-ERN.neg.fisher.txt.gz" &
bsub -J fisher -o fisher.out -e fisher.err -P $PROJECTID -K "python $BIN/fisher_test.py TS-ERP.pos.txt.gz TS-ERN.pos.txt.gz | gzip -c > TS-ERP_to_TS-ERN.pos.fisher.txt.gz" &
wait

# convert fisher test results to bed format
pksites=$HOME/projects/polya/results/common/polya_sites/PK.sites.c13.bed.gz
bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS.pos.fisher.txt.gz $pksites | bedtools sort -i - > NBT_to_TS.pos.fisher.bed" &
bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS.neg.fisher.txt.gz $pksites | bedtools sort -i - > NBT_to_TS.neg.fisher.bed" &
bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS-ERP.neg.fisher.txt.gz $pksites | bedtools sort -i - > NBT_to_TS-ERP.neg.fisher.bed" &
bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS-ERP.pos.fisher.txt.gz $pksites | bedtools sort -i - > NBT_to_TS-ERP.pos.fisher.bed" &
bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS-ERN.pos.fisher.txt.gz $pksites | bedtools sort -i - > NBT_to_TS-ERN.pos.fisher.bed" &
bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py NBT_to_TS-ERN.neg.fisher.txt.gz $pksites | bedtools sort -i - > NBT_to_TS-ERN.neg.fisher.bed" &
bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py TS-ERP_to_TS-ERN.pos.fisher.txt.gz $pksites | bedtools sort -i - > TS-ERP_to_TS-ERN.pos.fisher.bed" &
bsub -J fisher2bed -o fisher2bed.out -e fisher2bed.err -P $PROJECTID -K "python $BIN/visualize_fisher_shifts.py TS-ERP_to_TS-ERN.neg.fisher.txt.gz $pksites | bedtools sort -i - > TS-ERP_to_TS-ERN.neg.fisher.bed" &
wait

# convert bed to bb
bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $CHROM_SIZES NBT_to_TS.pos.fisher.bed" &
bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $CHROM_SIZES NBT_to_TS.neg.fisher.bed" &
bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $CHROM_SIZES NBT_to_TS-ERP.neg.fisher.bed" &
bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $CHROM_SIZES NBT_to_TS-ERP.pos.fisher.bed" &
bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $CHROM_SIZES NBT_to_TS-ERN.pos.fisher.bed" &
bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $CHROM_SIZES NBT_to_TS-ERN.neg.fisher.bed" &
bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $CHROM_SIZES TS-ERP_to_TS-ERN.pos.fisher.bed" &
bsub -J bed2bb -o bed2bb.out -e bed2bb.err -P $PROJECTID -K "bed2bb.py --type bed12 $CHROM_SIZES TS-ERP_to_TS-ERN.neg.fisher.bed" &
wait

# copy bbs over to the hub
cp *.bb $HUB/hg19
