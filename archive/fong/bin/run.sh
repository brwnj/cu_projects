#! /usr/bin/env bash
#BSUB -J miso_driver
#BSUB -e miso_driver.%J.err
#BSUB -o miso_driver.%J.out
#BSUB -q normal

for bam in /vol2/home/fonga/Bentley/Illumina_data/*.bam; do
    python /vol1/home/brownj/projects/fong/bin/misofix.py \
        --bam $bam \
        --index /vol1/home/brownj/ref/hg19/miso_isoforms \
        --out /vol2/home/fonga/Bentley/$(basename $bam _L003_R1_001_filtered.ProperPairs_sorted.bam) \
        --queue_name normal
done