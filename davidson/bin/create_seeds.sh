#!/usr/bin/env bash
#BSUB -J find_seeds[1-6]
#BSUB -e find_seeds.%J.%I.err
#BSUB -o find_seeds.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>2] rusage[mem=2] span[hosts=1]"
#BSUB -n 1

<<DOC
create tags and find seeds among the reads.
DOC

set -o nounset -o pipefail -o errexit -x

sample=$LSB_JOBINDEX

bin=$HOME/devel/iSSAKE
len=35
data=$HOME/projects/davidson/data
tags=$data/trbv.tags.fa
reads=$data/20120924/$sample.trm.fa.gz
seeds=$data/20120924/$sample.$len.seeds.fa.gz

python $bin/create_tags.py -v -l $len $tags \
    | python $bin/find_seeds.py -v - $reads \
    | gzip -c > $seeds