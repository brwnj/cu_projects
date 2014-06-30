#! /usr/bin/env bash
#BSUB -J trim[1-46]
#BSUB -e trim.%J.%I.err
#BSUB -o trim.%J.%I.out
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -q normal
#BSUB -n 1
#BSUB -P hits-clip

<<DOC
trim peaks by sample using my script
DOC

set -o nounset -o errexit -o pipefail -x

# rows of 10
sampleids=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30
            MP31 MP34 MP35 MP36 MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG
            MP42.TCGA MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11 
            PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42
            PK51 PK52 PK53 PK54 helaa helab) # do not have peaks yet PK61 PK62 PK63)
sample=${sampleids[$(($LSB_JOBINDEX - 1))]}
bin=$HOME/projects/hits-clip/bin
results=$HOME/projects/hits-clip/results/common/samples/$sample
gd=$HOME/projects/hits-clip/data/combined.rmd.genomedata

python $bin/trim_peaks.py -v -w 50 $results/$sample.peaks.rmd.bed.gz $gd \
    | gzip -c \
    > $results/$sample.peaks.rmd.trimmed.bed.gz