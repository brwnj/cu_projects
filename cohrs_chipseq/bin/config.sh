#!/usr/bin/env bash
set -o nounset

<<DOC
H3K9acB1_S5 H3K9ac  1
H3KacE1_S7  H3K9ac  2
H3PANA1_S1  H3PAN   1
H3PAND1_S3  H3PAN   2
INPUT_S2    input   1
INPUT_S4    input   2
RIGGC1_S6   R IgG   1
RIGGF1_S8   R IgG   2
DOC

SAMPLES=(H3K9acB1_S5 H3KacE1_S7 H3PANA1_S1 H3PAND1_S3 INPUT_S2 INPUT_S4 RIGGC1_S6 RIGGF1_S8)
RESULTS=$HOME/projects/cohrs_chipseq/results/common
DATA=$HOME/projects/cohrs_chipseq/data/common
NOVOIDX=$HOME/ref/hg19/hg19.9606.novoidx
