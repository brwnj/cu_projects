#!/usr/bin/env bash
#BSUB -J mask_fastq[1-44]
#BSUB -e mask_fastq.%J.%I.err
#BSUB -o mask_fastq.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

<<DOC
run jay's mask_fastq.py on all hits samples'
DOC

set -o nounset -o errexit -o pipefail -x

samples=(idx0 MP1 MP10 MP11 MP2 MP20
            MP21 MP22 MP23 MP24 MP30
            MP31 MP34 MP35 MP36 MP38
            MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG
            MP42.TCGA MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA
            MP45.ACTG MP45.TCGA MP7 MP9 PK11
            PK12 PK21 PK22 PK23 PK24
            PK31 PK32 PK33 PK41 PK42
            PK51 PK52 PK53 PK54)
sample=${samples[$LSB_JOBINDEX]}

bam=$HOME/projects/hits-clip/results/common/samples/$sample/$sample.rmd.bam
fastq=$HOME/projects/hits-clip/data/common/$(echo $sample | cut -d. -f1)/$sample.fastq.gz
masked_fastq=$sample.masked.fastq.gz

pyscript=$HOME/devel/cliptools/cliptools/mask_fastq.py

python $pyscript -v --debug $bam $fastq | gzip -c > $masked_fastq