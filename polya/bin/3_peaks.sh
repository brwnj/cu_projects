#! /usr/bin/env bash
#BSUB -J peaks[1-27]
#BSUB -e peaks.%J.%I.err
#BSUB -o peaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P pillai_kabos_polya

<<DOC
call peaks using macs2.
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/polya/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
result=$RESULT/$sample
bam=$result/$sample.bam
out=$result/$sample

# the output from macs2
peak=${out}_peaks.bed
narrowpeak=${out}_peaks.narrowPeak
summit=${out}_summits.bed
xls=${out}_peaks.xls

# some peaks were extending outside of genomic coords
clipped_peak=${out}_peaks.bed.clipped
clipped_summit=${out}_summits.bed.clipped

if [[ ! -f $peak ]]; then
    macs2 callpeak -t $bam -n $out --keep-dup auto \
        --nomodel --bw 30 --shiftsize 5 --call-summits
    bedClip $peak $CHROM_SIZES $clipped_peak
    bedClip $summit $CHROM_SIZES $clipped_summit
    mv $clipped_peak $peak
    mv $clipped_summit $summit
    rm -f $xls
fi
