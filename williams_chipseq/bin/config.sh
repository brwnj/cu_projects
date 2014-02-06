#!/usr/bin/env bash
set -o nounset

SAMPLES=(
3B5_hela_1
3B5_hela_2
3B5_hela_3
3B5_mitotic_hela_1
3B5_mitotic_hela_2
3B5_mitotic_hela_3
hela_input
mitotic_hela_input
polII_hela_1
polII_hela_2
polII_hela_3
polII_mitotic_hela_1
polII_mitotic_hela_2
polII_mitotic_hela_3
SC184_hela_1
SC184_hela_2
SC184_hela_3
SC184_mitotic_hela_1
SC184_mitotic_hela_2
SC184_mitotic_hela_3
)

PI=williams_chipseq
GENOME=hg19
DATA=$HOME/projects/williams_chipseq/data/common
RESULTS=$HOME/projects/williams_chipseq/results/common
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
SIZES=$HOME/ref/hg19/hg19.sizes
HUB=$RESULTS/hub

if [[ ! -d $HUB ]]; then
    mkdir -p $HUB
fi
