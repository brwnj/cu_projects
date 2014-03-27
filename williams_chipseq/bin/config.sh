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
3B5_hela
3B5_mitotic_hela
polII_hela
polII_mitotic_hela
SC184_hela
SC184_mitotic_hela
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

declare -A COLORS
COLORS=([3B5_hela_1]="0,109,44"
        [3B5_hela_2]="35,139,69"
        [3B5_hela_3]="65,174,118"
        [3B5_hela]="0,68,27"
        [3B5_mitotic_hela_1]="8,81,156"
        [3B5_mitotic_hela_2]="33,113,181"
        [3B5_mitotic_hela_3]="66,146,198"
        [3B5_mitotic_hela]="8,48,107"
        [polII_hela_1]="231,41,138"
        [polII_hela_2]="206,18,86"
        [polII_hela_3]="152,0,67"
        [polII_hela]="103,0,31"
        [polII_mitotic_hela_1]="221,52,151"
        [polII_mitotic_hela_2]="174,1,126"
        [polII_mitotic_hela_3]="122,1,119"
        [polII_mitotic_hela]="73,0,106"
        [SC184_hela_1]="128,125,186"
        [SC184_hela_2]="106,81,163"
        [SC184_hela_3]="84,39,143"
        [SC184_hela]="63,0,125"
        [SC184_mitotic_hela_1]="239,59,44"
        [SC184_mitotic_hela_2]="203,24,29"
        [SC184_mitotic_hela_3]="165,15,21"
        [SC184_mitotic_hela]="103,0,13"
        [hela_input]="82,82,82"
        [mitotic_hela_input]="115,115,115")

declare -A UCSC_GROUP
UCSC_GROUP=([3B5_hela_1]=s3B5_hela
            [3B5_hela_2]=s3B5_hela
            [3B5_hela_3]=s3B5_hela
            [3B5_hela]=s3B5_hela
            [3B5_mitotic_hela_1]=s3B5_mitotic_hela
            [3B5_mitotic_hela_2]=s3B5_mitotic_hela
            [3B5_mitotic_hela_3]=s3B5_mitotic_hela
            [3B5_mitotic_hela]=s3B5_mitotic_hela
            [polII_hela_1]=polII_hela
            [polII_hela_2]=polII_hela
            [polII_hela_3]=polII_hela
            [polII_hela]=polII_hela
            [polII_mitotic_hela_1]=polII_mitotic_hela
            [polII_mitotic_hela_2]=polII_mitotic_hela
            [polII_mitotic_hela_3]=polII_mitotic_hela
            [polII_mitotic_hela]=polII_mitotic_hela
            [SC184_hela_1]=SC184_hela
            [SC184_hela_2]=SC184_hela
            [SC184_hela_3]=SC184_hela
            [SC184_hela]=SC184_hela
            [SC184_mitotic_hela_1]=SC184_mitotic_hela
            [SC184_mitotic_hela_2]=SC184_mitotic_hela
            [SC184_mitotic_hela_3]=SC184_mitotic_hela
            [SC184_mitotic_hela]=SC184_mitotic_hela
            [hela_input]=hela
            [mitotic_hela_input]=mitotic_hela)
