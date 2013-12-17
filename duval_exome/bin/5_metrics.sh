#!/usr/bin/env bash
#BSUB -J metrics[1-16]
#BSUB -e metrics.%J.%I.err
#BSUB -o metrics.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P duval_exome

<<DOC
collect alignment and targeted alignment metrics
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/duval_exome/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
bam=$RESULTS/$sample/$sample.realign.bam
coveragebed=$RESULTS/$sample/$sample.coverage.bed.gz

if [[ ! -f $coveragebed ]]; then
    bedtools coverage -abam $bam -b $TARGETS | gzip -c > $coveragebed
fi

# coverage
coveragestats=$RESULTS/$sample/$sample.coverage_stats.txt
if [[ ! -f $coveragestats ]]; then
    zcat $coveragebed | cut -f4 | stats > $coveragestats
fi

# fraction covered
fracstats=$RESULTS/$sample/$sample.frac_stats.txt
if [[ ! -f $fracstats ]]; then
    zcat $coveragebed | cut -f7 | stats > $fracstats
fi

# process those files in python to come up with a comprehensive table
