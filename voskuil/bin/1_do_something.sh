#!/usr/bin/env bash
#BSUB -J jobname[1-2]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 4

<<DOC
genomic coords of all TAs into bed

TA characterization
chr start stop ta_name feature position_in_gene read_count

1   1   2   ta_1    gene    32%    4

gene characterization
chr start stop gene_name number_of_TAs number_of_TAs_with_reads

$ wc -l ta.positions.bed 
   91021 ta.positions.bed
   
$ awk '$3=="gene"' bps.gff > bps.genes.gff

$ bedtools intersect -wa -a ta.positions.bed -b bps.genes.gff | wc -l
   71457

$ awk 'BEGIN{OFS=FS="\t"}{lengt=$5-$4-1;$4=int($4+(lengt*.05));$5=int($5-(lengt*.2));print}' bps.genes.gff > bps.genes.lethalinsert.gff

$ bedtools intersect -wa -a ta.positions.bed -b bps.genes.lethalinsert.gff | wc -l
   52170

91021 TAs genome-wide
71457 are located within genes
52170 are located after the first 5% and before the last 20% of the gene's length

DOC

set -o nounset -o pipefail -o errexit -x

SAMPLES=(idx0 sample1 sample2)
SAMPLE=${SAMPLES[$LSB_JOBINDEX]}