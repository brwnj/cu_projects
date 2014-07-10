#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit -x

<<DOC
I don't see any pacbio reads for sample:
AAA164N20.TAGCTT.fastq
DOC

threads=40
email=jmbrown@bigelow.org

src=/mnt/scgc/data/PacBio/src
illumina=/mnt/scgc/data/PacBio/Illumina_data
pacbio=/mnt/scgc/data/PacBio/PacBio_reads/batch2
output=/mnt/scgc/data/PacBio/coassemblies
options="--email $email --threads $threads"

python $src/coassemble_shuffled.py $options AAA076-E06 $illumina/AAA076E06.GCCAAT.fastq $pacbio/PB0771_filtered_subreads.fasta $output/AAA076-E06
python $src/coassemble_shuffled.py $options AAA160-J20 $illumina/AAA160J20.GCCAAT.fastq $pacbio/PB0772_filtered_subreads.fasta $output/AAA160-J20
python $src/coassemble_shuffled.py $options AAA164-A08 $illumina/AAA164A08.CTTGTA.fastq $pacbio/PB0773_filtered_subreads.fasta $output/AAA164-A08
python $src/coassemble_shuffled.py $options AAA168-E21 $illumina/AAA168E21.TAGCTT.fastq $pacbio/PB0774_filtered_subreads.fasta $output/AAA168-E21
python $src/coassemble_shuffled.py $options AAA164-B23 $illumina/AAA164B23.ACTTGA.fastq $pacbio/PB0833_filtered_subreads.fasta $output/AAA164-B23
python $src/coassemble_shuffled.py $options AAA164-I21 $illumina/AAA164I21.TGACCA.fastq $pacbio/PB0834_filtered_subreads.fasta $output/AAA164-I21
python $src/coassemble_shuffled.py $options AAA160-C11 $illumina/AAA160C11.GCCAAT.fastq $pacbio/PB0835_filtered_subreads.fasta $output/AAA160-C11
python $src/coassemble_shuffled.py $options AAA164-A21 $illumina/AAA164A21.ACTTGA.fastq $pacbio/PB0836_filtered_subreads.fasta $output/AAA164-A21

pacbio=/mnt/scgc/data/PacBio/PacBio_reads/batch1
python $src/coassemble_shuffled.py $options AAA164M04 $illumina/AAA164M04.CGATGT.fastq $pacbio/AAA164M04.trimmed.fasta $output/AAA164M04
python $src/coassemble_shuffled.py $options AAA164P11 $illumina/AAA164P11.GATCAG.fastq $pacbio/AAA164P11.trimmed.fasta $output/AAA164P11
