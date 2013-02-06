#!/usr/bin/env bash
#BSUB -J tcr_repertoire[1-6]
#BSUB -e tcr_repertoire.%J.%I.err
#BSUB -o tcr_repertoire.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>2] rusage[mem=2] span[hosts=1]"
#BSUB -n 1

<<DOC
create tags and find seeds among the reads.
DOC

set -o nounset -o pipefail -o errexit -x

sample=$LSB_JOBINDEX

data=$HOME/projects/davidson/data/20120924
fasta=$data/$sample.trm.fa.gz
seeds=$data/$sample.35.seeds.fa.gz

iSSAKE -f $fasta -s $seeds -m 30 -b $sample