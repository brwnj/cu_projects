#!/usr/bin/env bash
#BSUB -J reprocessing
#BSUB -e reprocessing.%J.err
#BSUB -o reprocessing.%J.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1
#BSUB -P hits-clip

<<DOC
The CLIP samples to filter are:

PK11, PK31 and PK15 (MCF-7)
PK21, PK41 and PK52 (BT474)
P24, PK42 and PK54 (MDA-231)

I think we should do this slightly differently than with the newer samples:

    I'd take each aligned read, extend it 25bp downstream (rather than 100bp),
    and look for the presence of any of the following sequences in the extended
    alignment:

        GTGTCA
        GTGCCA
        GTGTCT
        GTCTCA
        GTGACA

    Discard any reads with matches in the extended alignment, and then call
    peaks as before (including combined peaks for each cell line).
DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(PK11 PK31 PK51 PK21 PK41 PK52 PK24 PK42 PK54)

# sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

fasta=$HOME/ref/hg19/hg19.fa
fasta_dir=$HOME/ref/hg19/fa_per_chr
sizes=$HOME/ref/hg19/hg19.sizes
results=/vol1/home/brownj/projects/hits-clip/results/20140317

# bam1=/vol1/home/brownj/projects/hits-clip/results/common/$sample/$sample.bam
# bam2=/vol1/home/brownj/projects/hits-clip/results/common/$sample/$sample.rmd.bam
#
# outbam1=$results/${sample}_filtered.bam
# outbam2=$results/${sample}_filtered.rmd.bam
#
# python ~/projects/hits-clip/bin/scripts/filter_bam.py -m all -d GTGTCA -d GTGCCA -d GTGTCT -d GTCTCA -d GTGACA -b 25 $bam1 $outbam1 $fasta
# python ~/projects/hits-clip/bin/scripts/filter_bam.py -m all -d GTGTCA -d GTGCCA -d GTGTCT -d GTCTCA -d GTGACA -b 25 $bam2 $outbam2 $fasta

# build a new genomedata archive
bam2gd.py -p hits-clip $sizes $fasta_dir *.bam

# call peaks

# combine peaks across replicates
# add to some hub somewhere
