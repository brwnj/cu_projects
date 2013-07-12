#! /usr/bin/env bash
#BSUB -J conifer
#BSUB -o analyze.%J.out
#BSUB -e analyze.%J.err
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P shaikh
#BSUB -q normal

<<doc
run conifer on tamim's samples.
doc

set -o nounset -o errexit -o pipefail -x

results=$HOME/projects/shaikh/results/common
if [[ ! -d $results/analysis ]]; then
    mkdir -p $results/analysis
fi
if [[ ! -d $results/analysis/plots ]]; then
    mkdir -p $results/analysis/plots
fi
bin=$HOME/opt/conifer_v0.2.2

probes=$bin/probes.txt
scree_plot=$results/analysis/tamim_scree.png
hdf5=$results/analysis/tamim.hdf5
calls=$results/analysis/tamim_calls.txt

# create the hdf5
if [[ ! -f $hdf5 ]]; then
    # first you should plot screeplot to find svd cutoff
    python $bin/conifer.py analyze \
        --probes $probes \
        --rpkm_dir $results/rpkm \
        --plot_scree $scree_plot \
        --svd 5 \
        --output $hdf5 
fi

# do the testing
if [[ ! -f $calls ]]; then
    python $bin/conifer.py call \
        --input $hdf5 \
        --output $calls
fi

# create plots
python $bin/conifer.py plotcalls \
    --input $hdf5 \
    --calls $calls \
    --outputdir $results/analysis/plots
