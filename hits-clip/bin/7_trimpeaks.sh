#! /usr/bin/env bash
#BSUB -J trim_peaks[1-19]
#BSUB -e trim_peaks.%J.%I.err
#BSUB -o trim_peaks.%J.%I.out
#BSUB -q normal

<<DOC
of the combined peaks, find the max intensity of each peak and expand to the size
of the ago footprint.
DOC

group=${SAMPLE_GROUPS[$(($LSB_JOBINDEX - 1))]}
replicates=${GROUP_REPLICATES[$(($LSB_JOBINDEX - 1))]}
replicate_count=$(($(grep -o " " <<<$replicates | wc -l) + 1))
results=$RESULTS/$group

for strand in pos neg; do
    tracks=""
    untrimmed=$results/$group.$strand.peaks.bed.gz
    trimmed=$results/$group.$strand.trimmed.peaks.bed.gz
    for sample in $replicates; do
        tracks="$tracks -t $sample.rmd_${strand}"
    done
    if [[ ! -f $trimmed ]]; then
        python $SRC/trim_peaks.py $tracks $untrimmed $GENOMEDATA | bedtools sort -i - | gzip -c > $trimmed
    fi
done
