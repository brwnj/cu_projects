#! /usr/bin/env bash
#BSUB -J align[1-22]
#BSUB -e align.%J.%I.err
#BSUB -o align.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P pillai_kabos_polya

<<DOC
Trim the UMI from the FASTQ, align trimmed reads using Novoalign suppressing 
all reads that align more than once, then remove UMI duplicates from the
alignment.
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/polya/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

unprocessed_fastq=$DATA/$sample.fq.gz
fastq=$DATA/$sample.umi.fq.gz
# trim the UMI
if [[ ! -f $fastq ]]; then
    umitools process_fastq $unprocessed_fastq NNNNNV | gzip -c > $fastq
fi

results=$RESULT/$sample
if [[ ! -d $results ]]; then
    mkdir -p $results
fi

umibam=$RESULT/$sample/$sample.UMIs_not_removed.bam
bam=$results/$sample.bam

# align the reads
if [[ ! -f $umibam ]]; then
    novoalign -d $NOVOIDX -f $fastq -o SAM -n 50 -s 5 -l 20 -r None -c 4 -k \
        2> $sample.alignment.txt \
        | samtools view -ShuF4 - \
        | samtools sort -o - $sample.temp -m 9500000000 \
        > $umibam
    samtools index $umibam
    # create bw
    bam2bw $umibam $CHROM_SIZES $PROJECTID TRUE
fi
# process the UMIs
if [[ ! -f $bam ]]; then
    umitools process_bam $umibam $bam
    samtools index $bam
    # create bw
    bam2bw $bam $CHROM_SIZES $PROJECTID TRUE 
fi