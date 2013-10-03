#! /usr/bin/env bash
#BSUB -J counts[1-69]
#BSUB -e counts.%J.%I.err
#BSUB -o counts.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

<<DOC
get the counts. overwrite existing count files.
DOC

slop=2

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
sites=$RESULT/polya_sites/${sample:0:2}.sites.c13.bed.gz
slopsitesneg=$RESULT/polya_sites/${sample:0:2}.sites.c13.slop.$slop.neg.bed.gz
slopsitespos=$RESULT/polya_sites/${sample:0:2}.sites.c13.slop.$slop.pos.bed.gz
results=$RESULT/$sample

# this adds slop to the polyA sites
if [[ ! -f $slopsitesneg ]]; then
    bedtools slop -b $slop -i $sites -g $CHROM_SIZES \
        | awk '$6 == "+"' \
        | bedtools sort -i - \
        | gzip -c > $slopsitesneg
fi
if [[ ! -f $slopsitespos ]]; then
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
                print full[1]"\t"full[2]"\t"$0}' \
        | cut -f1,2,4 \
        | sort -k1,1 -k2,2n \
        | gzip -c > $countsout 
done

# old method

# for strand in pos neg; do
#     bedg=$results/$sample.$strand.bedgraph.gz
#     out=$results/$sample.$strand.counts.txt.gz
#     # if [[ ! -f $out ]]; then
#     python $BIN/read_counts.py $bedg $sites $CHROM_SIZES | gzip -c > $out
#     # fi
# done
