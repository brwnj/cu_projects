#! /usr/bin/env bash
#BSUB -J polyac[1-10]
#BSUB -e polyac.%J.%I.err
#BSUB -o polyac.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

<<DOC
run tests on defined comparisons
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/polya/bin/config.sh

# comparisona=(MP54 MP54 MP54 MP54 MP55 MP55 MP55 MP56 MP56 MP57)
# comparisonb=(MP55 MP56 MP57 MP58 MP56 MP57 MP58 MP57 MP58 MP58)
comparisona=(MP59 MP59 MP59 MP59 MP60 MP60 MP60 MP61 MP61 MP62)
comparisonb=(MP60 MP61 MP62 MP63 MP61 MP62 MP63 MP62 MP63 MP63)
samplea=${comparisona[$(($LSB_JOBINDEX - 1))]}
sampleb=${comparisonb[$(($LSB_JOBINDEX - 1))]}

for strand in pos neg; do
    bga=$RESULT/$samplea/$samplea.$strand.bedgraph.gz
    bgb=$RESULT/$sampleb/$sampleb.$strand.bedgraph.gz
    test_result=$RESULT/testing_results/${samplea}_${sampleb}.$strand.txt

    python $BIN/polyac.py $REFBED $bga $bgb > $test_result
done