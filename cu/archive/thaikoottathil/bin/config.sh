#!/usr/bin/env bash
set -o nounset -o

SAMPLES=(10 11 12 14 15 16 18 22 23 24 3 7)                                 #12
# 1 5 9 13 19 20 21 17
TAGS=(STAT6_2h STAT6_2h_c cJun_2h cJun_2h_c STAT6_24h STAT6_24h_c cJun_24h cJun_24h_c)
PROJECT=$HOME/projects/thaikoottathil
RESULTS=$PROJECT/results/common
DATA=$PROJECT/data
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
BACKGROUND=$RESULTS/IgG
