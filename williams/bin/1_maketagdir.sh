#!/usr/bin/env bash
#BSUB -J maketagdir
#BSUB -e maketagdir.%J.%I.err
#BSUB -o maketagdir.%J.%I.out
#BSUB -q normal
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams/bin/config.sh

H3K4Me3_1=$RESULTS/H3K4Me3_1
if [[ ! -d $H3K4Me3_1 ]]; then
    echo "makeTagDirectory $H3K4Me3_1 $RESULTS/1FE/1FE.bam" | bsez maketagdir -P $PI
fi

H3K4Me3_2=$RESULTS/H3K4Me3_2
if [[ ! -d $H3K4Me3_2 ]]; then
    echo "makeTagDirectory $H3K4Me3_2 $RESULTS/2FE/2FE.bam" | bsez maketagdir -P $PI
fi

H3K4Me3_3=$RESULTS/H3K4Me3_3
if [[ ! -d $H3K4Me3_3 ]]; then
    echo "makeTagDirectory $H3K4Me3_3 $RESULTS/4MNPM/4MNPM.bam" | bsez maketagdir -P $PI
fi

AP2=$RESULTS/AP2
if [[ ! -d $AP2 ]]; then
    echo "makeTagDirectory $AP2 $RESULTS/3FE/3FE.bam" | bsez maketagdir -P $PI
fi

INPUT=$RESULTS/INPUT
if [[ ! -d $INPUT ]]; then
    echo "makeTagDirectory $INPUT $RESULTS/5FE_Input/5FE_Input.bam" | bsez maketagdir -P $PI
fi
