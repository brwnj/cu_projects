#!/usr/bin/env bash
#BSUB -J maketagdir
#BSUB -e maketagdir.%J.%I.err
#BSUB -o maketagdir.%J.%I.out
#BSUB -q normal
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -P thaikoottathil

set -o nounset -o pipefail -o errexit -x
pi=thaikoottathil
source $HOME/projects/$pi/bin/config.sh

stat62h=$RESULTS/STAT6_2h
if [[ ! -d $stat62h ]]; then
    echo "makeTagDirectory $stat62h $RESULTS/7/7.bam $RESULTS/15/15.bam" | bsez maketagdir -P $pi
fi

stat62hc=$RESULTS/STAT6_2h_c
if [[ ! -d $stat62hc ]]; then
    echo "makeTagDirectory $stat62hc $RESULTS/3/3.bam $RESULTS/11/11.bam" | bsez maketagdir -P $pi
fi

cjun2h=$RESULTS/cJun_2h
if [[ ! -d $cjun2h ]]; then
    echo "makeTagDirectory $cjun2h $RESULTS/16/16.bam" | bsez maketagdir -P $pi
fi

cjun2hc=$RESULTS/cJun_2h_c
if [[ ! -d $cjun2hc ]]; then
    echo "makeTagDirectory $cjun2hc $RESULTS/12/12.bam" | bsez maketagdir -P $pi
fi

stat624h=$RESULTS/STAT6_24h
if [[ ! -d $stat624h ]]; then
    echo "makeTagDirectory $stat624h $RESULTS/18/18.bam" | bsez maketagdir -P $pi
fi

stat624hc=$RESULTS/STAT6_24h_c
if [[ ! -d $stat624hc ]]; then
    echo "makeTagDirectory $stat624hc $RESULTS/22/22.bam" | bsez maketagdir -P $pi
fi

cjun24h=$RESULTS/cJun_24h
if [[ ! -d $cjun24h ]]; then
    echo "makeTagDirectory $cjun24h $RESULTS/24/24.bam" | bsez maketagdir -P $pi
fi

cjun24hc=$RESULTS/cJun_24h_c
if [[ ! -d $cjun24hc ]]; then
    echo "makeTagDirectory $cjun24hc $RESULTS/23/23.bam" | bsez maketagdir -P $pi
fi

background=$RESULTS/IgG
if [[ ! -d $background ]]; then
    echo "makeTagDirectory $background $RESULTS/14/14.bam $RESULTS/10/10.bam" | bsez maketagdir -P $pi
fi
