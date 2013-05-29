#!/usr/bin/env bash
#BSUB -J quantify[1-20]
#BSUB -e quantify.%J.%I.err
#BSUB -o quantify.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P duval_canine

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/canine/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fq=$DATA/$sample.trm.fastq.gz
results=$RESULTS/$sample
out=$results/$sample.quant.txt.gz
if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $out ]]; then
    python ~/devel/iqseq/iqseq.py quantify -c 15 -m 3 $fq | gzip -c > $out
fi
