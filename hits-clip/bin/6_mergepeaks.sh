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

GROUPS=(
BMEC
BT474
BT474estr
BT474herc
HELA
hMSC
HS27A
HS5
HUVEC
MCF7
MCF7estr
MDA231
MCF7_E_72
MCF7_TAM_72
MCF7_TAMRES
)
GROUP_REPLICATES=(
"MP42.ACTG MP45.ACTG MP45.TCGA"
"PK21 PK41 PK52"
"PK22 PK53"
"PK23"
"helaa helab"
"MP36 MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA"
"MP1 MP21 MP35"
"MP2 MP20 MP34"
"MP24 MP38"
"PK11 PK31 PK51"
"PK12 PK32"
"PK24 PK42 PK54"
"PK61"
"PK62"
"PK63"
)

group=${GROUPS[$(($LSB_JOBINDEX - 1))]}
replicates=${GROUP_REPLICATES[$(($LSB_JOBINDEX - 1))]}
replicate_count=$(($(grep -o " " <<<$replicates | wc -l) + 1))

################################################
# PK62.rmd_pos.peaks.qv.passed_filter.bed.gz
################################################
results=$RESULTS/$group

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

for strand in pos neg; do
    bedfiles=""
    out=$results/$group.$strand.some.extension.gz
    if [[ ! -f $out ]]; then
        for sample in $replicates; do
            bedfiles="$bediles $RESULTS/$sample/$sample.rmd_${strand}.peaks.qv.passed_filter.bed.gz"
        done
        if [[ replicate_count == 1 ]]; then
            gunzip $bedfiles
            bedClip ${bedfiles/.gz} $SIZES $out
            gzip $out
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
