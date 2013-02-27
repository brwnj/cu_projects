#! /usr/bin/env bash
#BSUB -J "conifer[1-9,16-24]%4"
#BSUB -o conifer.%J.%I.out
#BSUB -e conifer.%J.%I.err
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -P shaikh

<<doc
run conifer on tamim's samples
trios only
doc

set -o nounset -o errexit -o pipefail -x

sample=$LSB_JOBINDEX
results=$HOME/projects/shaikh/results/common/rpkm
if [[ ! -d $results ]]; then
    mkdir -p $results
fi
bam=/vol1/home/gowank/Projects/TamimShaikh/bams/$sample.NewRG.bam
bin=$HOME/opt/conifer_v0.2.2
rpkm=$results/$sample.rpkm.txt

python $bin/conifer.py rpkm \
    --probes $bin/probes.txt \
    --input $bam \
    --output $rpkm