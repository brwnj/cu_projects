#!/usr/bin/env bash

reads=$HOME/projects/duval_exome/data/20130819
results=$HOME/projects/duval_exome/data/common

echo "cat $reads/1025_TGACCA_L002_R1_001.fastq.gz $reads/1025_TGACCA_L003_R1_001.fastq.gz > $results/1025_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/1025_TGACCA_L002_R2_001.fastq.gz $reads/1025_TGACCA_L003_R2_001.fastq.gz > $results/1025_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/113_CTTGTA_L006_R1_001.fastq.gz $reads/113_CTTGTA_L007_R1_001.fastq.gz > $results/113_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/113_CTTGTA_L006_R2_001.fastq.gz $reads/113_CTTGTA_L007_R2_001.fastq.gz > $results/113_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/166408_AAAGCA_L006_R1_001.fastq.gz $reads/166408_AAAGCA_L007_R1_001.fastq.gz > $results/166408_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/166408_AAAGCA_L006_R2_001.fastq.gz $reads/166408_AAAGCA_L007_R2_001.fastq.gz > $results/166408_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/172160_CAAAAG_L006_R1_001.fastq.gz $reads/172160_CAAAAG_L007_R1_001.fastq.gz > $results/172160_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/172160_CAAAAG_L006_R2_001.fastq.gz $reads/172160_CAAAAG_L007_R2_001.fastq.gz > $results/172160_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/22_GATCAG_L006_R1_001.fastq.gz $reads/22_GATCAG_L007_R1_001.fastq.gz > $results/22_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/22_GATCAG_L006_R2_001.fastq.gz $reads/22_GATCAG_L007_R2_001.fastq.gz > $results/22_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/353_GAAACC_L006_R1_001.fastq.gz $reads/353_GAAACC_L007_R1_001.fastq.gz > $results/353GAAACC_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/353_GAAACC_L006_R2_001.fastq.gz $reads/353_GAAACC_L007_R2_001.fastq.gz > $results/353GAAACC_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/353_GGCTAC_L006_R1_001.fastq.gz $reads/353_GGCTAC_L007_R1_001.fastq.gz > $results/353GGCTAC_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/353_GGCTAC_L006_R2_001.fastq.gz $reads/353_GGCTAC_L007_R2_001.fastq.gz > $results/353GGCTAC_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/36_CAGATC_L002_R1_001.fastq.gz $reads/36_CAGATC_L003_R1_001.fastq.gz > $results/36_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/36_CAGATC_L002_R2_001.fastq.gz $reads/36_CAGATC_L003_R2_001.fastq.gz > $results/36_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/400_TTAGGC_L002_R1_001.fastq.gz $reads/400_TTAGGC_L003_R1_001.fastq.gz > $results/400_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/400_TTAGGC_L002_R2_001.fastq.gz $reads/400_TTAGGC_L003_R2_001.fastq.gz > $results/400_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/500_TAGCTT_L006_R1_001.fastq.gz $reads/500_TAGCTT_L007_R1_001.fastq.gz > $results/500_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/500_TAGCTT_L006_R2_001.fastq.gz $reads/500_TAGCTT_L007_R2_001.fastq.gz > $results/500_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/522_AAACAT_L006_R1_001.fastq.gz $reads/522_AAACAT_L007_R1_001.fastq.gz > $results/522_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/522_AAACAT_L006_R2_001.fastq.gz $reads/522_AAACAT_L007_R2_001.fastq.gz > $results/522_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/868_ACTTGA_L002_R1_001.fastq.gz $reads/868_ACTTGA_L003_R1_001.fastq.gz > $results/868_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/868_ACTTGA_L002_R2_001.fastq.gz $reads/868_ACTTGA_L003_R2_001.fastq.gz > $results/868_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/Bililey_CGATGT_L002_R1_001.fastq.gz $reads/Bililey_CGATGT_L003_R1_001.fastq.gz > $results/Bililey_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/Bililey_CGATGT_L002_R2_001.fastq.gz $reads/Bililey_CGATGT_L003_R2_001.fastq.gz > $results/Bililey_R2.fastq.gz" | bsez concat -P duval
echo "cat $reads/KTCC_ATCACG_L002_R1_001.fastq.gz $reads/KTCC_ATCACG_L003_R1_001.fastq.gz > $results/KTCC_R1.fastq.gz" | bsez concat -P duval
echo "cat $reads/KTCC_ATCACG_L002_R2_001.fastq.gz $reads/KTCC_ATCACG_L003_R2_001.fastq.gz > $results/KTCC_R2.fastq.gz" | bsez concat -P duval

echo "zcat $reads/730_GCCAAT_L002_R1_001.fastq.gz $reads/730_GCCAAT_L003_R1_001.fastq.gz | gzip -c > $results/730_R1.fastq.gz" | bsez concat -P duval
echo "zcat $reads/730_GCCAAT_L002_R2_001.fastq.gz $reads/730_GCCAAT_L003_R2_001.fastq.gz | gzip -c > $results/730_R2.fastq.gz" | bsez concat -P duval
echo "zcat $reads/735_ACAGTG_L002_R1_001.fastq.gz $reads/735_ACAGTG_L003_R1_001.fastq.gz | gzip -c > $results/735_R1.fastq.gz" | bsez concat -P duval
echo "zcat $reads/735_ACAGTG_L002_R2_001.fastq.gz $reads/735_ACAGTG_L003_R2_001.fastq.gz | gzip -c > $results/735_R2.fastq.gz" | bsez concat -P duval
