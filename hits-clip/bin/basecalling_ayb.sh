#! /usr/bin/env bash

<<DOC
run ayb with varying number of model iterations
DOC

set -o nounset -o pipefail -o errexit -x

NOMASK=R57
MASK=R6I14C36
ITERATIONS=5
# L4T1101-2316

for (( i = 1101; i < 2317; i++ )); do
    RUNSCRIPT=ayb.${i}.sh
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J ayb.control.${i}" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1] select[mem>8] rusage[mem=8]\"" >> $RUNSCRIPT
    echo "#BSUB -e %J.err" >> $RUNSCRIPT
    echo "#BSUB -q night" >> $RUNSCRIPT
    echo "

AYB --dataformat cif --format fastq --blockstring $NOMASK --output ayb.${ITERATIONS} --niter ${ITERATIONS} --loglevel debug --input ~/projects/hits-clip/data/intensity_files --runfolder L1T${i} s_+" >> $RUNSCRIPT

   bsub < $RUNSCRIPT
done