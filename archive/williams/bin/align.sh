#!/usr/bin/env bash
#BSUB -J novoalign
#BSUB -e novoalign.%J.err
#BSUB -o novoalign.%J.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 4

<<DOC
align to cs indx incorporating latest dbsnp info
DOC

set -o nounset -o pipefail -o errexit -x

novoidx=$HOME/ref/mm9/mm9.colorspace.10090.novoidx
csfasta=$HOME/projects/williams/data/20130108/all_trimmed_F3.csfasta.gz
results=$HOME/projects/williams/results/common/all
bam=$results/all.bam

novoalignCS -d $novoidx -F CSFASTAnQV -f $csfasta -o SAM -r A 20 \
    -c 4 -k -s 5 -l 25 2> align_stats.txt \
    | samtools view -ShuF4 - \
    | samtools sort -o - all.temp -m 9500000000 > $bam