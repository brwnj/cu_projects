#!/usr/bin/env bash
#BSUB -J align[1-28]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 12
#BSUB -P koch

set -o nounset -o pipefail -o errexit -x

SAMPLES=(
Douglas_A1_ATCACG_L004_R1_001
Douglas_A2_CGATGT_L004_R1_001
Douglas_A3_TTAGGC_L004_R1_001
Douglas_B1_TGACCA_L004_R1_001
Douglas_B2_ACAGTG_L004_R1_001
Douglas_B3_GCCAAT_L004_R1_001
Douglas_C1_CAGATC_L008_R1_001
Douglas_C2_ACTTGA_L008_R1_001
Douglas_C3_GATCAG_L008_R1_001
Douglas_D1_TAGCTT_L008_R1_001
Douglas_D2_GGCTAC_L008_R1_001
Douglas_D3_CTTGTA_L008_R1_001
MUT2_6plus_34plus_TGACCA_L001_R1_001
MUT2_6plus_CGATGT_L001_R1_001
MUT3_6plus_34plus_AGTCAA_L001_R1_001
MUT3_6plus_GTGAAA_L001_R1_001
WT2_6plus_34plus_ATGTCA_L001_R1_001
WT2_6plus_CTTGTA_L001_R1_001
WT3_6plus_34plus_CCGTCC_L001_R1_001
WT3_6plus_GTCCGC_L001_R1_001
1_AAGGGA_L003_R1_001
2_CCTTCA_L003_R1_001
3_GGACCC_L003_R1_001
4_TTCAGC_L003_R1_001
5_AAGACG_L003_R1_001
6_CCTCGG_L003_R1_001
7_GGATGT_L003_R1_001
8_TTCGCT_L003_R1_001
)

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
fastq=/vol1/home/brownj/projects/edwards_mapping/data/20140327/$sample.trimmed.fastq.gz
results=/vol1/home/brownj/projects/edwards_mapping/results/common/$sample
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

bam=$results/${sample}.bam
if [[ ! -f $bam ]]; then
    STAR --runThreadN 12 \
        --genomeDir ~analysiscore/genomes/reference/mm10/mm10_GencodeM2 \
        --readFilesIn $fastq \
        --readFilesCommand zcat \
        --outFileNamePrefix $results/${sample}_ \
        --outFilterMultimapNmax 2 \
        --sjdbGTFfile ~analysiscore/genomes/reference/mm10/mm10_GencodeM2/gencode.vM2.annotation.gtf

    # clean up the STAR output
    sam=$results/${sample}_Aligned.out.sam

    samtools view -ShuF4 $sam | samtools sort -@ 12 -m 16G - ${bam/.bam}
    samtools index $bam
    rm $sam
fi
