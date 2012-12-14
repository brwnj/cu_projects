#! /usr/bin/env bash
#BSUB -J solid2fastq[1-]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

samples=(idx0
DK13_Dunn_P18_F3
DK13_Dunn_P18_F5
DK15_Dunn_Tumor2_F3
DK15_Dunn_Tumor2_F5
DK16_Dunn_Tumor3_F3
DK16_Dunn_Tumor3_F5
DK17_Dunn_Tumor4_F3
DK17_Dunn_Tumor4_F5
DK19_DLM8_Tumor1_F3
DK19_DLM8_Tumor1_F5
DK1_DLM8_F3
DK1_DLM8_F5
DK20_DLM8_Tumor2_F3
DK20_DLM8_Tumor2_F5
DK21_DLM8_Tumor3_F3
DK21_DLM8_Tumor3_F5
DK23_Dunn_P9_F3
DK23_Dunn_P9_F5
DK25_Dunn_P10_F3
DK25_Dunn_P10_F5
DK27_DLM8_F3
DK27_DLM8_F5
DK29_DLM8_2DG_F3
DK29_DLM8_2DG_F5
DK2_DLM8_F3
DK2_DLM8_F5
DK31_DLM8_2DG_F3
DK31_DLM8_2DG_F5
DK6_DLM8_2DG_F3
DK6_DLM8_2DG_F5
Dunn13_F3
Dunn13_F5
Dunn1_F3
Dunn1_F5
Dunn2_F3
Dunn2_F5
Dunn3_F3
Dunn3_F5
)

solid2fastq -o prefix test.csfasta test.qual 