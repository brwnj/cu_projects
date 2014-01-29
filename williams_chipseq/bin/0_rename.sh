#!/usr/bin/env bash
#BSUB -J rename
#BSUB -e rename.%J.err
#BSUB -o rename.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

<<DOC

DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

raw=$HOME/projects/williams_chipseq/data/20140124/140114_7001413_0114_AC3EA2ACXX
data=$HOME/projects/williams_chipseq/data/common

mv $raw/10_CACACA_L001_R1_001.fastq.gz $data/hela_input.fastq.gz
mv $raw/11_TTCAGC_L001_R1_001.fastq.gz $data/polII_mitotic_hela_1.fastq.gz
mv $raw/12_CACCTC_L008_R1_001.fastq.gz $data/3B5_mitotic_hela_1.fastq.gz
mv $raw/13_AAGGGA_L008_R1_001.fastq.gz $data/SC184_mitotic_hela_1.fastq.gz
mv $raw/14_AAGACG_L001_R1_001.fastq.gz $data/polII_mitotic_hela_2.fastq.gz
mv $raw/15_GTGGCC_L008_R1_001.fastq.gz $data/3B5_mitotic_hela_2.fastq.gz
mv $raw/16_CCTTCA_L008_R1_001.fastq.gz $data/SC184_mitotic_hela_2.fastq.gz
mv $raw/17_CCTCGG_L001_R1_001.fastq.gz $data/polII_mitotic_hela_3.fastq.gz
mv $raw/18_TGTTGC_L008_R1_001.fastq.gz $data/3B5_mitotic_hela_3.fastq.gz
mv $raw/19_GGACCC_L008_R1_001.fastq.gz $data/SC184_mitotic_hela_3.fastq.gz
mv $raw/1_AAGGGA_L001_R1_001.fastq.gz $data/polII_hela_1.fastq.gz
mv $raw/20_TTCAGC_L008_R1_001.fastq.gz $data/mitotic_hela_input.fastq.gz
mv $raw/2_GTGTTA_L008_R1_001.fastq.gz $data/3B5_hela_1.fastq.gz
mv $raw/3_GGATGT_L001_R1_001.fastq.gz $data/SC184_hela_1.fastq.gz
mv $raw/4_CCTTCA_L001_R1_001.fastq.gz $data/polII_hela_2.fastq.gz
mv $raw/5_TGTGAA_L008_R1_001.fastq.gz $data/3B5_hela_2.fastq.gz
mv $raw/6_TTCGCT_L001_R1_001.fastq.gz $data/SC184_hela2.fastq.gz
mv $raw/7_GGACCC_L001_R1_001.fastq.gz $data/polII_hela_3.fastq.gz
mv $raw/8_ACAAAC_L008_R1_001.fastq.gz $data/3B5_hela_3.fastq.gz
mv $raw/9_ACACGA_L001_R1_001.fastq.gz $data/SC184_hela_3.fastq.gz
