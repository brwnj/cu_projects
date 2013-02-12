#!/usr/bin/env bash
#BSUB -J process_assemblies[1-6]
#BSUB -e process_assemblies.%J.%I.err
#BSUB -o process_assemblies.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>32] rusage[mem=32] span[hosts=1]"
#BSUB -n 1

<<DOC
create tags and find seeds among the reads.
1, 3, and 5 are alpha
2, 4, and 6 are beta
DOC

set -o nounset -o pipefail -o errexit -x

sample=$LSB_JOBINDEX

data=$HOME/projects/davidson/data/20120924
results=$HOME/projects/davidson/results/20130206
bin=$HOME/devel/iSSAKE

rem=$(( $sample % 2 ))
if [ $rem -eq 0 ]; then
  # even number
  target=$data/trbj.fa
else
  # odd number
  target=$data/traj.fa
fi

contigs=$results/$sample.contigs
annotated_fasta=$results/$sample.filtered.fa
metadata=$results/$sample.metadata
aligned_fasta=$results/$sample.aligned.fa

exonerate \
    -q $contigs \
    -t $target \
    --bestn 1 \
    --ryo ">%qi|%ti\n%qs" \
    --showalignment FALSE \
    --showvulgar FALSE \
    | grep -v "Command line:\|Hostname:\|-- completed" \
    > $annotated_fasta
python $bin/reads2meta.py $annotated_fasta > $metadata
muscle -in $annotated_fasta -out $aligned_fasta

gzip $contigs $annotated_fasta $metadata $aligned_fasta
