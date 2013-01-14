#! /usr/bin/env bash
#BSUB -J trim.umi[1-12]
#BSUB -e trim.umi%J.%I.err
#BSUB -o trim.umi%J.%I.out
#BSUB -q short

<<DOC
trim the umi from the read and incorporate it into the read name.
DOC

set -o nounset -o pipefail -o errexit -x

# SAMPLES=(idx0 MP51 MP52 MP53 PK61 PK62)
# SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}

samples=(idx0 100-3 86-7 93-1 93-3 m+c m+e2
            nbt29 nbt39 nbt89 ts21 ts28 ts57)
sample=${samples[${LSB_JOBINDEX}]}

bin=$HOME/devel/umitools/umitools
reads=$HOME/projects/polya/data/20130114/${sample}.fq.gz
trim_result=$HOME/projects/polya/data/20130114/${sample}.umi.fq.gz

python $bin/process_fastq.py $reads NNNNNV | gzip -c > $trim_result