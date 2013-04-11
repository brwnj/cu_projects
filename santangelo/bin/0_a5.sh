#! /usr/bin/env bash
#BSUB -J a5
#BSUB -e a5.%J.err
#BSUB -o a5.%J.out
#BSUB -q bigmem
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -P santangelo

<<DOC
assemble bacterial genome using A5.
DOC

set -o nounset -o pipefail -o errexit -x
# source $HOME/projects/santangelo/bin/config.sh

# samples=(KW128)
# sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

sample=KW128
reads=$HOME/projects/santangelo/data
r1=$reads/${sample}_R1.fastq
r2=$reads/${sample}_R2.fastq
results=$HOME/projects/santangelo/results/common/$sample
out=$results/$sample
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

~/opt/ngopt_20130326/bin/a5_pipeline.pl $r1 $r2 $sample