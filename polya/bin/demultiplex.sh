#!/usr/bin/env bash
#BSUB -J demultiplex[1-2]
#BSUB -e demultiplex.%J.%I.err
#BSUB -o demultiplex.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 6

<<DOC
demultiplex using ea-tools' fastq-multx
DOC

set -o nounset -o errexit -x

SAMPLES=(idx0 6 7)
SAMPLE=${SAMPLES[$LSB_JOBINDEX]}

# had to trim 1 base from the 3' end due to AYB cycles
# then complement the bases of the given indexes

nbt29 TGTCAC
ts21 CGGTTA
nbt39 GTCTAG
ts28 TGAACT
nbt89 CTAGTC
ts57 ATCGAA

fastq-multx -B L006_complement.txt -m 2 6.fastq.gz -o %.fq