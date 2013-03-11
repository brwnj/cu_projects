#! /usr/bin/env bash
#BSUB -J genomedata
#BSUB -e genomedata.%J.err
#BSUB -o genomedata.%J.out
#BSUB -P hits-clip

set -o nounset -o errexit -o pipefail -x

<<DOC
Create or add tracks to genomedata archive.

To add tracks:
    genomedata-open-data
    genomedata-load-data
    genomedata-close-data
DOC

# will include duplicate and non-duplicate bams
archive=$HOME/projects/hits-clip/data/combined.rmd.genomedata
# sizes=$HOME/ref/hg18/hg18.sizes
# fasta=$HOME/ref/hg18/fa_per_chr
# bams="$HOME/projects/hits-clip/results/common/samples/*/*.bam"
# 
# bam2gd -o $archive $sizes $fasta $bams

# this should be able to add tracks as well

# It's worth noting that I changed the way files were named in later iterations
# of this pipeline.

# tracks to add
samples="PK61 PK62 PK63"
results=$HOME/projects/hits-clip/results/common/samples
for sample in $samples; do
    for strand in pos neg; do
        # note the order of strand and rmd extension between bedgraph and track
        bedgraph=$results/$sample/$sample.rmd.$strand.bedgraph.gz
        track=$sample.$strand.rmd
        genomedata-open-data $archive $track
        zcat $bedgraph | genomedata-load-data $archive $track
    done
done
genomedata-close-data $archive