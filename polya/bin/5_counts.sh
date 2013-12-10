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

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
slopsitesneg=$POLYASITES/${sample:0:2}.test_sites.slop.$SLOP.neg.bed.gz
slopsitespos=$POLYASITES/${sample:0:2}.test_sites.slop.$SLOP.pos.bed.gz
results=$RESULTS/$sample

for strand in pos neg; do

    bedg=$results/$sample.$strand.bedgraph.gz
    countsout=$results/$sample.$strand.counts.txt.gz
    slopsites=$POLYASITES/${sample:0:2}.test_sites.slop.$SLOP.$strand.bed.gz

    python $BIN/read_counts.py $bedg $slopsites | gzip -c > $countsout

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
