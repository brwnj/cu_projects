#! /usr/bin/env bash
#### BSUB -J ayb[1-2]
#### BSUB -e ayb.%J.%I.err
#### BSUB -o ayb.%J.%I.out
#### BSUB -q normal
#### BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#### BSUB -n 12

<<DOC
call bases using AYB
DOC

set -o nounset -o errexit -x

# lanes=("idx0" "6" "7")
# lane=${lanes[${LSB_JOBINDEX}]}

bs=R6I14C37
cifs=/vol1/home/brownj/projects/polya/data/20121210
out=/vol1/home/brownj/projects/polya/results/20121214

# AYB -b $bs -d cif -l debug -o $out -p 12 -i $cifs -r L${lane}T1101-2316 --format fastq

# you need a for loop here because some seg fault
for (( i = 1101; i < 2320; i++ )); do
    RUNSCRIPT=ayb.${i}.sh
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J ayb.${i}" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1] select[mem>8] rusage[mem=8]\"" >> $RUNSCRIPT
    echo "#BSUB -e ayb.%J.err" >> $RUNSCRIPT
    echo "#BSUB -o ayb.%J.out" >> $RUNSCRIPT
    echo "#BSUB -n 10" >> $RUNSCRIPT
    echo "#BSUB -q normal" >> $RUNSCRIPT
    echo "

AYB -b $bs -d cif -l debug -o $out -p 10 -i $cifs -r L7T${i} --format fastq s_+
" >> $RUNSCRIPT

    bsub < $RUNSCRIPT
done