#! /usr/bin/env bash
#BSUB -J trim.umi[1-5]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q short

<<DOC
trim the umi from the read and incorporate it into the read name.
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(idx0 MP51 MP52 MP53 PK61 PK62)
SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}

BIN=$HOME/projects/hits-clip/bin
DATA=$HOME/projects/hits-clip/results/20120717/ayb.mask.5

python $BIN/umitools-trim-umi.py $DATA/$SAMPLE.fq.gz NNNNNV | gzip -c > $SAMPLE.umi.fq.gz