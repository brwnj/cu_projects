#! /usr/bin/env bash
#BSUB -J genomedata
#BSUB -e genomedata.%J.err
#BSUB -o genomedata.%J.out

set -o nounset -o errexit -o pipefail -x

# will include duplicate and non-duplicate bams
archive=$HOME/projects/hits-clip/data/gd_20130227
sizes=$HOME/ref/hg18/hg18.sizes
fasta=$HOME/ref/hg18/fa_per_chr
bams="$HOME/projects/hits-clip/results/common/samples/*/*.bam"

bam2gd -o $archive $sizes $fasta $bams