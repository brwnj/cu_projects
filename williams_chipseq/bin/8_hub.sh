#!/usr/bin/env bash
#BSUB -J hub[1-20]
#BSUB -e hub.%J.%I.err
#BSUB -o hub.%J.%I.out
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
posbg=$RESULTS/$sample/${sample}_pos.bedgraph
negbg=$RESULTS/$sample/${sample}_neg.bedgraph
posbw=$RESULTS/$sample/${sample}_pos.bw
negbw=$RESULTS/$sample/${sample}_neg.bw

if [[ ! -f $posbg ]]; then
    makeUCSCfile $tagdir -strand + -norm $normalize > $posbg
    bedGraphToBigWig $posbg $SIZES $posbw
    makeUCSCfile $tagdir -strand - -norm $normalize > $negbg
    bedGraphToBigWig $negbg $SIZES $negbw
fi


# create bigbeds of the peaks
awk -t '$1!~/^#/{print $2,$3,$4,$1,$6,$5}' ../common/3B5_hela_1/peaks.txt | bedtools sort -i - | head
