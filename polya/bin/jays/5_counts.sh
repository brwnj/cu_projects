#! /usr/bin/env bash
#BSUB -J counts[1-6]
#BSUB -e counts.%J.%I.err
#BSUB -o counts.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

# XXX: should make this an option
slop=2

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sites=$RESULT/polya_sites/${sample:0:2}.sites.c13.bed.gz
slopsitesneg=$RESULT/polya_sites/${sample:0:2}.sites.c13.slop.$slop.neg.bed.gz
slopsitespos=$RESULT/polya_sites/${sample:0:2}.sites.c13.slop.$slop.pos.bed.gz
results=$RESULT/$sample

# XXX neg strand sites need to be compared with positive strand data.
# write out sites flipped from their strand sign

if [[ ! -f $slopsitesneg ]]; then
    # this adds slop to the polyA sites
    bedtools slop -b $slop -i $sites -g $CHROM_SIZES \
        | awk '$6 == "+"' \
        | bedtools sort -i - \
        | gzip -c > $slopsitesneg
    bedtools slop -b $slop -i $sites -g $CHROM_SIZES \
        | awk '$6 == "-"' \
        | bedtools sort -i - \
        | gzip -c > $slopsitespos
fi

for strand in pos neg; do

    bedg=$results/$sample.$strand.bedgraph.gz
    countsout=$results/$sample.$strand.counts.txt.gz
    slopsites=$RESULT/polya_sites/${sample:0:2}.sites.c13.slop.$slop.$strand.bed.gz

    # this mess writes out counts in gene/site/count format
    bedtools map -c 4 -o max -null 0 -a $slopsites -b $bedg \
        | awk '{split($4, symbol, "|"); split(symbol[1], gene, ".");
                print gene[3]":"$4"\t"$7}' \
        | awk '{split($1, full, ":"); split(full[2], site, ".");
                print full[1]"\t"full[2]"\t"$0}'  \
        | cut -f1,2,4 \
        | sort -k1,1 -k2,2n \
        | gzip -c > $countsout 
done
