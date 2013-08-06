#!/usr/bin/env bash
set -o nounset

<<DOC
1, 2, 4 histone marker H3K4Me3
3 transcription factor AP-2
5 input
DOC

PI=williams_chipseq
SAMPLES=(1FE 2FE 3FE 4MNPM 5FE_Input)
TAGS=(H3K4Me3_1 H3K4Me3_2 AP2 H3K4Me3_3)
STYLES=(histone histone factor histone)
DATA=$HOME/projects/williams/data/20130806
RESULTS=$HOME/projects/williams/results/common
NOVOIDX=$HOME/ref/mm10/mm10.nix
BACKGROUND=$RESULTS/INPUT
