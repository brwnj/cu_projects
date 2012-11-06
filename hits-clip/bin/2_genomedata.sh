#! /usr/bin/env bash
#BSUB -J loaddata[1-25]
#BSUB -e %J.err
#BSUB -o %J.out

# job array to parallelize by chr
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
        chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18
        chr19 chr20 chr21 chr22 chrX chrY chrM)
CHROMIDX=$(expr $LSB_JOBINDEX - 1)
CHROM=${CHROMS[$CHROMIDX]}

# writes the arrays
# SAMPLEIDS="MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP34 MP35 MP36 MP7 MP9 
#             PK11 PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 
#             PK53 PK54 MP30 MP31 MP38 MP40 MP41 MP39.ACTG MP39.TCGA MP42.ACTG 
#             MP42.TCGA MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG 
#             MP45.TCGA"
SAMPLEIDS="helaa helab"
STRANDS="pos neg"
DATA="/vol1/home/brownj/projects/hits-clip/results/common/samples"

VARS=vars.$CHROM.sh
# required by load_chr.sh
echo "SEQDIR=/vol3/home/jhessel/projects/encode/data/hg18/fasta" > $VARS
echo "GENOMEDATADIR=$CHROM.gd" >> $VARS
echo "BED_TRACKNAMES=($(for SAMPLE in $SAMPLEIDS; do for STRAND in $STRANDS; do echo -n "$SAMPLE.$STRAND "; done; done))" >> $VARS
echo "BED_FILENAMES=($(for SAMPLE in $SAMPLEIDS; do for STRAND in $STRANDS; do echo -n "$DATA/$SAMPLE/$SAMPLE.$STRAND.bedgraph.gz "; done; done))" >> $VARS

JOBNAME="genomedata_load.$CHROM"
# ~/opt/bin/load_chr_sequence.sh
bsub -J $JOBNAME -o %J.out -e %J.err "load_chr_sequence.sh $VARS $CHROM"