#! /usr/bin/env bash
#BSUB -J tophat_colorspace[1-19]
#BSUB -e tophat.%J.%I.err
#BSUB -o tophat.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 4

set -o nounset -o errexit -x

samples=(idx0
DK13_Dunn_P18
DK15_Dunn_Tumor2
DK16_Dunn_Tumor3
DK17_Dunn_Tumor4
DK19_DLM8_Tumor1
DK1_DLM8
DK20_DLM8_Tumor2
DK21_DLM8_Tumor3
DK23_Dunn_P9
DK25_Dunn_P10
DK27_DLM8
DK29_DLM8_2DG
DK2_DLM8
DK31_DLM8_2DG
DK6_DLM8_2DG
Dunn13
Dunn1
Dunn2
Dunn3
)

sample=${samples[$LSB_JOBINDEX]}

out=/vol1/home/brownj/projects/kamstock/results/common/$sample
bidx=/vol1/home/brownj/ref/mm9/mm9c
gtf=/vol1/home/brownj/ref/mm9/mm9.ncbi37.gtf
input=/vol1/home/brownj/projects/kamstock/data/20121203

# tophat2 --color --quals --read-realign-edit-dist 0 --bowtie1 -o $out \
#     -r 200 --mate-std-dev 50 -a 5 -p 4 --no-discordant --no-mixed \
#     --library-type fr-secondstrand -G $gtf -T $bidx \
#     $input/${sample}_F3_trim_F3.csfasta $input/${sample}_F5_trim_F5.csfasta \
#     $input/${sample}_F3_trim_F3_QV.qual $input/${sample}_F5_trim_F5_QV.qual

tophat2 --color --quals --read-realign-edit-dist 2 --bowtie1 -o $out \
    -r 200 --mate-std-dev 50 -a 5 -p 4 \
    --library-type fr-secondstrand -G $gtf -T $bidx \
    $input/${sample}_F3_trim_F3.csfasta $input/${sample}_F5_trim_F5.csfasta \
    $input/${sample}_F3_trim_F3_QV.qual $input/${sample}_F5_trim_F5_QV.qual