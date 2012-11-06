#! /usr/bin/env bash
<<DOC
Concatenate reads into one file.
DOC

set -o nounset -o pipefail -o errexit -x

echo "zcat 1_MD4xML_ACAGTG_L005_R1_001.fastq.gz 1_MD4xML_ACAGTG_L007_R1_001.fastq.gz | gzip -c > 1_MD4xML_ACAGTG.fastq.gz" | bsub-ez concat
echo "zcat 2_MD4B220_CAGATC_L005_R1_001.fastq.gz 2_MD4B220_CAGATC_L007_R1_001.fastq.gz | gzip -c > 2_MD4B220_CAGATC.fastq.gz" | bsub-ez concat
echo "zcat 3_ARS_A1_ATGTCA_L005_R1_001.fastq.gz 3_ARS_A1_ATGTCA_L007_R1_001.fastq.gz | gzip -c > 3_ARS_A1_ATGTCA.fastq.gz" | bsub-ez concat
echo "zcat 4_BL6B220_CGATGT_L005_R1_001.fastq.gz 4_BL6B220_CGATGT_L007_R1_001.fastq.gz | gzip -c > 4_BL6B220_CGATGT.fastq.gz" | bsub-ez concat
echo "zcat 5_BL6Alpha_GTGAAA_L005_R1_001.fastq.gz 5_BL6Alpha_GTGAAA_L007_R1_001.fastq.gz | gzip -c > 5_BL6Alpha_GTGAAA.fastq.gz" | bsub-ez concat
echo "zcat 6_BL6Gamma_GTCCGC_L005_R1_001.fastq.gz 6_BL6Gamma_GTCCGC_L007_R1_001.fastq.gz | gzip -c > 6_BL6Gamma_GTCCGC.fastq.gz" | bsub-ez concat
echo "zcat 7_BL6-FO_CCGTCC_L005_R1_001.fastq.gz 7_BL6-FO_CCGTCC_L007_R1_001.fastq.gz | gzip -c > 7_BL6-FO_CCGTCC.fastq.gz" | bsub-ez concat
echo "zcat 8_BL6-ACT_TGACCA_L005_R1_001.fastq.gz 8_BL6-ACT_TGACCA_L007_R1_001.fastq.gz | gzip -c > 8_BL6-ACT_TGACCA.fastq.gz" | bsub-ez concat