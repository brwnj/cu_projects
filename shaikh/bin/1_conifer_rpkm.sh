#! /usr/bin/env bash
#BSUB -J "conifer[1-43]%8"
#BSUB -o conifer.%J.%I.out
#BSUB -e conifer.%J.%I.err
#BSUB -R "span[hosts=1]"
#BSUB -n 1
#BSUB -P shaikh

<<DOC
run conifer on tamim's trio samples with taylor's families as added background.
DOC

set -o nounset -o errexit -o pipefail -x

# 10 samples per row
samples=(1 2 3 4 5 6 7 8 9 16
17 18 19 20 21 22 23 24 MT48_25 MT48_26
MT48_27 MT48_33 MT48_34 MT48_37 MT48_30 MT48_31 MT48_32 MT48_28 MT48_29 MT48_35
MT48_1 MT48_2 MT48_3 MT48_15 MT48_19 MT48_20 MT48_21 MT48_12 MT48_22 MT16_3
MT48_9 MT48_10 MT48_11)
sample=${samples[$(($LSB_JOBINDEX - 1))]}

if [[ $sample == MT* ]]; then
    data=/vol1/home/gowank/Projects/MattTaylor/16and48/WholeBAMs
    ext=for64.bam
else
    data=/vol1/home/gowank/Projects/TamimShaikh/bams
    ext=NewRG.bam
fi

bam=$data/$sample.$ext

results=$HOME/projects/shaikh/results/common/rpkm
bin=$HOME/opt/conifer_v0.2.2
rpkm=$results/$sample.rpkm.txt

if [[ ! -d $results ]]; then
    mkdir -p $results
fi
if [[ ! -f $rpkm ]]; then
    python $bin/conifer.py rpkm \
        --probes $bin/hg19_ens_gene.txt \
        --input $bam \
        --output $rpkm
fi
