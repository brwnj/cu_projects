#!/usr/bin/env bash
#BSUB -J merge_peaks[1-15]
#BSUB -e merge_peaks.%J.%I.err
#BSUB -o merge_peaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P hits-clip

<<DOC
Merges the replicate groups without regard to strand and strand specific. It
also filters the bed where start is greater than stop and clips coordinates
down to chromosome sizes.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/hits-clip/bin/config.sh

group=${SAMPLE_GROUPS[$((LSB_JOBINDEX - 1))]}
replicates=${GROUP_REPLICATES[$((LSB_JOBINDEX - 1))]}
replicate_count=$(($(grep -c -o " " <<<$replicates | wc -l) + 1))

results=$RESULTS/$group

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

for strand in pos neg; do
    bedfiles=""
    out=$results/$group.$strand.peaks.bed.gz
    if [[ ! -f $out ]]; then
        for sample in $replicates; do
            bedfiles="$bedfiles $RESULTS/$sample/$sample.rmd_${strand}.peaks.qv.passed_filter.bed.gz"
        done
        if [[ $replicate_count == 1 ]]; then
            gunzip $bedfiles
            bedClip ${bedfiles/.gz} $SIZES ${out/.gz}
            gzip -f ${out/.gz}
            gzip ${bedfiles/.gz}
        else
            toclip=${out/.gz/.clipme}
            peaktools-combine-replicates --verbose $bedfiles > $toclip
            bedClip $toclip $SIZES $out
            gzip $out
            rm -f $toclip
        fi
    fi
done
