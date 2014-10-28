#! /usr/bin/env bash
#BSUB -J cuffdiff[1-5]
#BSUB -e cuffdiff.%J.%I.err
#BSUB -o cuffdiff.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 4

set -o nounset -o errexit -x

<<DOC
Dunn vs Dunn Tumor
------------------------------------------------------------------------------
Dunn13,Dunn1,Dunn2,Dunn3 DK15_Dunn_Tumor2,DK16_Dunn_Tumor3,DK17_Dunn_Tumor4

Dunn vs DLM8
------------------------------------------------------------------------------
Dunn13,Dunn1,Dunn2,Dunn3 DK1_DLM8,DK2_DLM8,DK27_DLM8

Dlm8 vs DLM8 tumor
------------------------------------------------------------------------------
DK1_DLM8,DK2_DLM8,DK27_DLM8 DK19_DLM8_Tumor1,DK20_DLM8_Tumor2,DK21_DLM8_Tumor3

Dunn tumor vs DLM8 tumor
------------------------------------------------------------------------------
DK15_Dunn_Tumor2,DK16_Dunn_Tumor3,DK17_Dunn_Tumor4 DK19_DLM8_Tumor1,DK20_DLM8_Tumor2,DK21_DLM8_Tumor3

DLM8 vs. DLM8 2-DG
------------------------------------------------------------------------------
DK1_DLM8,DK2_DLM8,DK27_DLM8 DK29_DLM8_2DG,DK31_DLM8_2DG,DK6_DLM8_2DG
DOC

input=$HOME/projects/kamstock/results/common
comparisons=(idx0 "Dunn_vs_DunnTumor" "Dunn_vs_DLM8" "DLM8_vs_DLM8Tumor" "DunnTumor_vs_DLM8Tumor" "DLM8_vs_DLM82DG")
comparison=${comparisons[${LSB_JOBINDEX}]}

Dunn="$input/Dunn13/accepted_hits.bam,$input/Dunn1/accepted_hits.bam,$input/Dunn2/accepted_hits.bam,$input/Dunn3/accepted_hits.bam"
DunnTumor="$input/DK15_Dunn_Tumor2/accepted_hits.bam,$input/DK16_Dunn_Tumor3/accepted_hits.bam,$input/DK17_Dunn_Tumor4/accepted_hits.bam"
DLM8="$input/DK1_DLM8/accepted_hits.bam,$input/DK2_DLM8/accepted_hits.bam,$input/DK27_DLM8/accepted_hits.bam"
DLM8Tumor="$input/DK19_DLM8_Tumor1/accepted_hits.bam,$input/DK20_DLM8_Tumor2/accepted_hits.bam,$input/DK21_DLM8_Tumor3/accepted_hits.bam"
DLM82DG="$input/DK29_DLM8_2DG/accepted_hits.bam,$input/DK31_DLM8_2DG/accepted_hits.bam,$input/DK6_DLM8_2DG/accepted_hits.bam"

leftsamples=(idx0 $Dunn $Dunn $DLM8 $DunnTumor $DLM8)
rightsamples=(idx0 $DunnTumor $DLM8 $DLM8Tumor $DLM8Tumor $DLM82DG)

leftsample=${leftsamples[${LSB_JOBINDEX}]}
rightsample=${rightsamples[${LSB_JOBINDEX}]}

leftlabel=$(echo $comparison | cut -f1 -d_)
rightlabel=$(echo $comparison | cut -f3 -d_)
out=$HOME/projects/kamstock/results/common/$comparison
gtf=$HOME/ref/mm9/mm9.ncbi37.gtf
fasta=$HOME/ref/mm9/mm9.fa

cuffdiff -o $out -b $fasta -u -L $leftlabel,$rightlabel -p 4 \
    --library-type fr-secondstrand $gtf $leftsample $rightsample