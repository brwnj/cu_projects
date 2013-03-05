#!/usr/bin/env bash
#BSUB -J ayb
#BSUB -e ayb.%J.err
#BSUB -o ayb.%J.out
#BSUB -q normal
#BSUB -P pillai_kabos_polya

<<DOC
call bases using AYB
DOC

set -o nounset -o errexit -x

# lanes=("idx0" "6" "7")
# lane=${lanes[${LSB_JOBINDEX}]}
lane=3

# bs=R6I14C36 # barcode is 3'; UMI is 5'
bs=R6I14C36
cifs=$HOME/projects/polya/data/20130305
out=$HOME/projects/polya/results/20130305

# AYB -b $bs -d cif -l debug -o $out -p 12 -i $cifs -r L${lane}T1101-2316 --format fastq

# you need a for loop here because some seg fault
for (( i = 1101; i < 2320; i++ )); do
    RUNSCRIPT=ayb.${i}.$lane.sh
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J ayb.${i}.$lane" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1] select[mem>8] rusage[mem=8]\"" >> $RUNSCRIPT
    echo "#BSUB -e ayb.%J.$lane.err" >> $RUNSCRIPT
    echo "#BSUB -o ayb.%J.$lane.out" >> $RUNSCRIPT
    echo "#BSUB -n 10" >> $RUNSCRIPT
    echo "#BSUB -q normal" >> $RUNSCRIPT
    echo "#BSUB -P pillai_kabos_polya" >> $RUNSCRIPT
    echo "

AYB -b $bs -d cif -l debug -o $out -p 10 -i $cifs -r L${lane}T${i} --format fastq s+
" >> $RUNSCRIPT

    bsub < $RUNSCRIPT
done