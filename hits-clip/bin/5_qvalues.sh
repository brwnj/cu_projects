#!/usr/bin/env bash
#BSUB -J qvalues
#BSUB -e qvalues.%J.err
#BSUB -o qvalues.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_hitsclip

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/hits-clip/bin/config.sh

cutoff=.001

# more serial processing
for track in `genomedata-info tracknames_continuous $GENOMEDATA`; do
    sample=$(echo ${track/_*} | cut -f1 -d.)
    real=$RESULTS/$sample/$track.peaks.bed.gz
    null=$RESULTS/$sample/$track.shuffle.peaks.bed.gz
    all=$RESULTS/$sample/$track.peaks.qv.all.bed.gz
    filtered=$RESULTS/$sample/$track.peaks.qv.passed_filter.bed.gz
    summary=$RESULTS/$sample/$track.qvalue_summary.txt
    if [[ ! -f $filtered ]]; then
        peaktools-qvalues -v $real $null | gzip -c > $all
        awk -tc header -v CO=${cutoff} '$7<CO' $all | gzip -c > $filtered
        peaktools-qvalue-summary -t $cutoff $all > $summary
    fi
done
