#!/bin/sh
#NH - number of reported alignments that contains the query in the current record

#for sample in A_BC01 J_BC02 L_BC03 M_BC04: do qsub -v dir=$sample bam_process.sh

samtools view -F 0x04 -h /mnt/storage1/rnaseq/WT_Asthma_AJKL_Paired_End/$dir/output/pairing/wt.pe.bam
    | awk 'BEGIN {FS=OFS="\t"}
        (!/^@/){
            minus=and($2, 0x10);
            split($12, a, ":");
            $12=a[1]":"a[2]":"a[3]+1;
            print $0"\t"XS:A:(minus ? "-":"+")
            }
        (/^@/){
            print
        }' 
    | samtools view -bhS - 
    > /mnt/storage1/rnaseq/WT_Asthma_AJKL_Paired_End/$dir/output/pairing/wt.pe.xs.nh.bam
