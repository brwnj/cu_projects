#!/usr/bin/env bash
#BSUB -J process_assemblies[1-6]
#BSUB -e process_assemblies.%J.%I.err
#BSUB -o process_assemblies.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 1

<<DOC
create tags and find seeds among the reads.
1, 3, and 5 are alpha
2, 4, and 6 are beta
DOC

set -o nounset -o pipefail -o errexit -x

sample=$LSB_JOBINDEX

data=
results=
bin=

rem=$(( $sample % 2 ))
if [ $rem -eq 0 ]; then
  # even number
  target=trbj.fa
else
  # odd number
  target=traj.fa
fi

contigs=
exonerate_out=
annotated_fasta=
metadata=
aligned_fasta=

exonerate -q $contigs -t $target --bestn 1 --ryo ">%qi|%ti\n%qs" --showalignment FALSE --showvulgar FALSE > $exonerate_out
grep -v "Command line:\|Hostname:\|-- completed" $exonerate_out > $annotated_fasta
python $bin reads2meta.py $annotated_fasta > $metadata
muscle -in $annotated_fasta -out $aligned_fasta