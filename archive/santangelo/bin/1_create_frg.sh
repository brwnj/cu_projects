#! /usr/bin/env bash
#BSUB -J fastqToCA[1-7]
#BSUB -e fastqToCA.%J.%I.err
#BSUB -o fastqToCA.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>1] rusage[mem=1] span[hosts=1]"
#BSUB -n 1
#BSUB -P santangelo

<<DOC
Create FRG file for each sample.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/santangelo/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
r1=$READS/${sample}_R1.fastq
r2=$READS/${sample}_R2.fastq
frg=$READS/$sample.frg

fastqToCA -insertsize 800 150 -libraryname $sample -type sanger -mates $r1,$r2 > $frg