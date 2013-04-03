#! /usr/bin/env bash
#BSUB -J join_reads[1-6]
#BSUB -e join_reads.%J.%I.err
#BSUB -o join_reads.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

<<DOC
Join R1 and R2 into continuous sequence using SeqPrep.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
r1=$READS/${sample}_R1.filtered.trimmed.fastq.gz
r2=$READS/${sample}_R2.filtered.trimmed.fastq.gz
results=$RESULTS/$sample
r1out=$results/${sample}_R1.joined.fastq.gz
r2out=$results/${sample}_R2.joined.fastq.gz
joined=$results/$sample.joined.fastq.gz

if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $joined ]]; then
    SeqPrep -f $r1 -r $r2 -1 $r1out -2 $r2out -s $joined
fi