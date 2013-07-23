#!/usr/bin/env bash
#BSUB -J annotatepeaks[1-4]
#BSUB -e annotatepeaks.%J.%I.err
#BSUB -o annotatepeaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P thaikoottathil

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/thaikoottathil/bin/config.sh

as=(STAT6_2h cJun_2h STAT6_24h cJun_24h)
bs=(STAT6_2h_c cJun_2h_c STAT6_24h_c cJun_24h_c)
a=${as[$(($LSB_JOBINDEX - 1))]}
b=${bs[$(($LSB_JOBINDEX - 1))]}
peaks=$RESULTS/$a/peaks.txt
diffpeaks=$RESULTS/$a/peaks.diff.txt
adir=$RESULTS/$a
bdir=$RESULTS/$b
# all.motif is a concatenation of known and discovered motifs
motif=$RESULTS/${a}_motif/all.motif
out=$RESULTS/$a/peaks.annotated.txt
diffout=$RESULTS/$a/peaks.diff.annotated.txt

if [[ ! -f $motif ]]; then
    cat $RESULTS/${a}_motif/homerResults/motif*.motif > $motif
    cat $RESULTS/${a}_motif/knownResults/known*.motif >> $motif
fi
if [[ ! -f $out ]]; then
    # do we really need to annotate the peaks of the controls?
    annotatePeaks.pl $peaks hg19 -size 400 -d $adir $bdir -m $motif > $out
fi
if [[ ! -f $diffout ]]; then
    annotatePeaks.pl $diffpeaks hg19 -size 400 -d $adir $bdir -m $motif > $diffout
fi
