#!/usr/bin/env bash
#BSUB -J hubprep[1-20]
#BSUB -e hubprep.%J.%I.err
#BSUB -o hubprep.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

<<DOC
convert tag dirs to bigwigs
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
tagdir=$RESULTS/$sample
normalize=1.5e7


# create the bigwigs
posbg=$tagdir/${sample}_pos.bedgraph
negbg=$tagdir/${sample}_neg.bedgraph
posbw=$tagdir/${sample}_pos.bw
negbw=$tagdir/${sample}_neg.bw

if [[ ! -f $posbg ]]; then
    makeUCSCfile $tagdir -strand + -norm $normalize > $posbg
    bedGraphToBigWig $posbg $SIZES $posbw
    makeUCSCfile $tagdir -strand - -norm $normalize > $negbg
    bedGraphToBigWig $negbg $SIZES $negbw
fi


# create bigbeds of the peaks
peaks=$tagdir/peaks.txt
peakbed=$tagdir/${sample}_peaks.bed
peakbigbed=$tagdir/${sample}_peaks.bb

if [[ ! -f $peakbigbed ]]; then
    awk -t '$1!~/^#/{print $2,$3,$4,$1,"1",$5}' $peaks | bedtools sort -i - > $peakbed
    bedToBigBed -type=bed6 $peakbed $SIZES $peakbigbed
    gzip -f $peakbed
fi
