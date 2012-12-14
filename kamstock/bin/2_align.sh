#! /usr/bin/env bash
#BSUB -J solid2fastq[1-]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

mouse

solid2fastq -o sample_name sample.csfasta sample.qual

tophat2 + bowtie1?

