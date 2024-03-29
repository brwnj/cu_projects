#!/usr/bin/env bash
#BSUB -J tcr_repertoire[1,3,5]
#BSUB -e tcr_repertoire.%J.%I.err
#BSUB -o tcr_repertoire.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>12] rusage[mem=12] span[hosts=1]"
#BSUB -n 1
#BSUB -P davidson

<<DOC
create tags and find seeds among the reads.
1, 3, and 5 are alpha
2, 4, and 6 are beta
DOC

set -o nounset -o pipefail -o errexit -x

sample=$LSB_JOBINDEX

data=$HOME/projects/davidson/data/20120924
fasta=$data/$sample.trm.fa.gz
# seeds=$data/$sample.35.seeds.fa.gz

results=$HOME/projects/davidson/results/common
seeds=$results/$sample/$sample.contigs.gz

iSSAKE -f $fasta -s $seeds -b $sample
