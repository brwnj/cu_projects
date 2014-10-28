#!/usr/bin/env bash
set -o nounset

PI=spillman
SAMPLES=(E21PR1_minus E21PR1_plus E21PR2_minus E21PR2_plus E21PR3_minus #5
            E21PR3_plus E2Input_minus E2input_plus)                     #8
DATA=$HOME/projects/spillman/data/20130625
RESULTS=$HOME/projects/spillman/results/common
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
