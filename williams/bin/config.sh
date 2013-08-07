#!/usr/bin/env bash
set -o nounset

<<DOC
1, 2, 4 histone marker H3K4Me3
3 transcription factor AP-2
5 input

	                1FE	        2FE	        3FE	        4MNPM	    5FE_Input
Read Sequences	    26737040	27140651	27666285	33157316	9559600
Aligned	            26287725	26151794	20857409	32468329	9423154
Unique Alignment	24676204	23926175	16937105	30137603	7097320
Gapped Alignment	484871	    490475	    521705	    604628	    271115
Quality Filter	    15452	    13162	    6143	    18728	    2216
Homopolymer Filter	4644	    4579	    3575	    4841	    1045
DOC

PI=williams_chipseq
SAMPLES=(1FE 2FE 3FE 4MNPM 5FE_Input)
TAGS=(H3K4Me3_1 H3K4Me3_2 AP2 H3K4Me3_3)
STYLES=(histone histone factor histone)
DATA=$HOME/projects/williams/data/20130806
RESULTS=$HOME/projects/williams/results/common
NOVOIDX=$HOME/ref/mm10/mm10.nix
BACKGROUND=$RESULTS/INPUT
