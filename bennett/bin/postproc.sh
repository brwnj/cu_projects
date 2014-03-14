#!/usr/bin/env bash
#BSUB -J postproc[1-72]
#BSUB -e postproc.%J.%I.err
#BSUB -o postproc.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>6] rusage[mem=6] span[hosts=1]"
#BSUB -n 1
#BSUB -P bennett

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/bennett/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
bin=$HOME/projects/bennett/bin
results=$RESULTS/$sample

if [[ ! -d $results ]]; then
    mkdir -p $results
fi

imgt_aa=$results/${sample}_imgtaa.txt.gz
dists_result=$results/${sample}_unique_aa.tsv
stats=$results/${sample}_dists.tsv
plot=$results/${sample}_dist.pdf

# distributions
python $bin/aa_dists.py -m 6 -s .70 $imgt_aa 2> $stats > $dists_result
# plotting
Rscript $bin/plot_dist.R $dists_result $sample $plot
