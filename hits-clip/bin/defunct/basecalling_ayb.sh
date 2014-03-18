#! /usr/bin/env bash

<<DOC
run ayb with varying number of model iterations
DOC

set -o nounset -o errexit -x

MASK=R6I14C36
ITERATIONS=5
# L4T1101-2316
lowerbound=
upperbound=
cifs=$HOME/projects/polya/data/20121210/fsfsdfsdf

for (( i = $lowerbound; i < $upperbound; i++ )); do
    RUNSCRIPT=ayb.${i}.sh
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J ayb.${i}" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1] select[mem>8] rusage[mem=8]\"" >> $RUNSCRIPT
    echo "#BSUB -e %J.err" >> $RUNSCRIPT
    echo "#BSUB -q normal" >> $RUNSCRIPT
    echo "

AYB --dataformat cif --format fastq --blockstring $MASK --output ayb \
--niter $ITERATIONS --loglevel debug --input $cifs \
--runfolder L1T${i} s_+" >> $RUNSCRIPT

   bsub < $RUNSCRIPT
done