#!/usr/bin/env bash
#BSUB -J combinepeaks
#BSUB -e combine.%J.err
#BSUB -o combine.%J.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

<<DOC
combine replicate peaks
DOC

set -o nounset -o pipefail -o errexit -x

data=$HOME/projects/walter/results/common
bin=$HOME/devel/peaktools/peaktools

python $bin/combine_replicates.py -m 4 $data/*Inf/*tb2hg19.qvalues.gz | gzip -c > $data/inf.peaks.bed.gz
python $bin/combine_replicates.py -m 4 $data/*Uninf/*tb2hg19.qvalues.gz | gzip -c > $data/uninf.peaks.bed.gz