#! /usr/bin/env bash
#BSUB -J ayb[1-2]
#BSUB -e ayb.%J.%I.err
#BSUB -o ayb.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 12

<<DOC
call bases using AYB
DOC

set -o nounset -o errexit -x

lanes=("idx0" "6" "7")
lane=${lanes[${LSB_JOBINDEX}]}

bs=R6I14C37
cifs=/vol1/home/brownj/projects/polya/data/20121210
out=/vol1/home/brownj/projects/polya/results/20121214

AYB -b $bs -d cif -l debug -o $out -p 12 -i $cifs -r L${lane}T1101-2316 --format fastq