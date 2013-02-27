#! /usr/bin/env bash
#BSUB -J conifer
#BSUB -o analyze.%J.out
#BSUB -e analyze.%J.err
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P shaikh

<<doc
run conifer on tamim's samples
trios only
doc

set -o nounset -o errexit -o pipefail -x

results=$HOME/projects/shaikh/results/common
if [[ ! -d $results/analysis ]]; then
    mkdir -p $results/analysis
fi
bin=$HOME/opt/conifer_v0.2.2

python $bin/conifer.py analyze \
    --probes $bin/probes.txt \
    --rpkm_dir $results/rpkm \
    --output $results/analysis/trios.hdf5 \
    --write_svals $results/analysis/trios_singular_values.txt \
    --write_sd $results/analysis/trios_sd_values.txt


import matplotlib.pyplot as plt
import pylab as P
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle