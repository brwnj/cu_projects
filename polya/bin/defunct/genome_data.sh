#! /usr/bin/env bash
#BSUB -J load_data[1-25]
#BSUB -e load_data.%J.%I.err
#BSUB -o load_data.%J.%I.out

<<DOC
create genomedata archive with sequence data
DOC

set -o nounset -o pipefail -o errexit -x

# job array to parallelize by chr
CHROMS=(idx0 chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18
        chr19 chr20 chr21 chr22 chrX chrY chrM)
CHROM=${CHROMS[$LSB_JOBINDEX]}

# writes the arrays
SAMPLEIDS="MP51 MP52 MP53 PK61 PK62"
STRANDS="pos neg"
ANNOTATIONS="umi noumi"
EXTS="umi.bedgraph.gz bedgraph.gz"
DATA=/vol1/home/brownj/projects/polya/results/common

VARS=vars.$CHROM.sh
# required by load_chr.sh
echo "SEQDIR=/vol3/home/jhessel/projects/encode/data/hg18/fasta" > $VARS
echo "GENOMEDATADIR=$CHROM.gd" >> $VARS
echo "BED_TRACKNAMES=($(for SAMPLE in $SAMPLEIDS; do for ANNOTATION in $ANNOTATIONS; do echo -n "$SAMPLE.$ANNOTATION "; done; done))" >> $VARS
echo "BED_FILENAMES=($(for SAMPLE in $SAMPLEIDS; do for EXT in $EXTS; do echo -n "$DATA/$SAMPLE/$SAMPLE.$EXT "; done; done))" >> $VARS

JOBNAME="genomedata_load.$CHROM"
# ~/opt/bin/load_chr_sequence.sh
bsub -J $JOBNAME -o %J.out -e %J.err "load_chr_sequence.sh $VARS $CHROM"