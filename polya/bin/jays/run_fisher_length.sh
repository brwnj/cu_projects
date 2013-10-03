#! /usr/bin/env bash
#BSUB -J fisher.length
#BSUB -e fisher.length.%J.err
#BSUB -o fisher.length.%J.out
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

# XXX orig
# results="$RESULT/fisher_results"

results="$RESULT/fisher_results/flipped.test"
if [[ ! -d $results ]]; then
    mkdir $results
fi

for strand in pos neg; do

    shorta="$RESULT/PC3/PC3.$strand.counts.txt.gz"
    shortb="$RESULT/PC4/PC4.$strand.counts.txt.gz"

    # test flipped samples
    longb="$RESULT/PC5/PC5.$strand.counts.txt.gz"
    longa="$RESULT/PC6/PC6.$strand.counts.txt.gz"

    # XXX original way
    # longa="$RESULT/PC5/PC5.$strand.counts.txt.gz"
    # longb="$RESULT/PC6/PC6.$strand.counts.txt.gz"

    out="$results/fisher.length.$strand.tab.gz"
    python $BIN/fisher_test_length.py \
        --longa $longa \
        --longb $longb \
        --shorta $shorta \
        --shortb $shortb \
        --verbose  \
        | gzip -c > $out
done 

