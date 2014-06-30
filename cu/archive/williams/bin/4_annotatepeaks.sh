#!/usr/bin/env bash
#BSUB -J annotatepeaks[1-4]
#BSUB -e annotatepeaks.%J.%I.err
#BSUB -o annotatepeaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams/bin/config.sh

tag=${TAGS[$(($LSB_JOBINDEX - 1))]}

tagdir=$RESULTS/$tag
peaks=$tagdir/regions.txt
out=$RESULTS/$tag/regions.annotated.txt
if [[ $tag == "AP2" ]]; then
    peaks=$tagdir/peaks.txt
    out=$RESULTS/$tag/peaks.annotated.txt
fi

# all.motif is a concatenation of known and discovered motifs
motif=$RESULTS/${tag}_motif/all.motif

if [[ ! -f $motif ]]; then
    cat $RESULTS/${tag}_motif/homerResults/motif*.motif > $motif
    cat $RESULTS/${tag}_motif/knownResults/known*.motif >> $motif
fi
if [[ ! -f $out ]]; then
    annotatePeaks.pl $peaks mm10 -d $tagdir -m $motif > $out
fi
