#!/usr/bin/env bash
#BSUB -J 3prime_seed_extension
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 4
#BSUB -cwd /vol1/home/brownj/projects/davidson/results/common/1

<<DOC
run ssake given trav and trbv seeds
DOC

set -o nounset -o pipefail -o errexit -x

reads=/vol1/home/brownj/projects/davidson/data/20120924
out=/vol1/home/brownj/projects/davidson/results/common/1

seed_fa=$reads/tr_ab_v.fa
gzip_fasta=$reads/1.jnd.fa.gz
fasta=$out/1.jnd.fa

if [[ ! -f $fasta ]]; then
    zcat $gzip_fasta > $fasta
fi

SSAKE -f $fasta -s $seed_fa -m 40 -o 50 -r 0.8 -b s1 -p 1 -v 1 -d 200 -e 0.75 -k 10 -a 0.5 -x 50