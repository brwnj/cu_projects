#!/usr/bin/env bash
set -o nounset -o pipefail -o errexit -x

<<DOC
I don't see any pacbio reads for sample:
AAA164N20.TAGCTT.fastq
DOC

threads=30
email=jmbrown@bigelow.org

src=/mnt/scgc/data/PacBio/src

illumina=/mnt/scgc/data/PacBio/Assemblies
pacbio=/mnt/scgc/data/PacBio/PacBio_reads/batch2
output=/mnt/scgc/data/PacBio/coassemblies
options="--email $email --threads $threads"

# pacbio coassemblies with kmernorm and complexity filtering
# python $src/coassemble_shuffled.py $options AAA076-E06 $illumina/AAA076E06.GCCAAT/AAA076E06.GCCAAT.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/PB0771_filtered_subreads.fasta $output/AAA076-E06
# python $src/coassemble_shuffled.py $options AAA160-J20 $illumina/AAA160J20.GCCAAT/AAA160J20.GCCAAT.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/PB0772_filtered_subreads.fasta $output/AAA160-J20
# python $src/coassemble_shuffled.py $options AAA164-A08 $illumina/AAA164A08.CTTGTA/AAA164A08.CTTGTA.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/PB0773_filtered_subreads.fasta $output/AAA164-A08
# python $src/coassemble_shuffled.py $options AAA168-E21 $illumina/AAA168E21.TAGCTT/AAA168E21.TAGCTT.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/PB0774_filtered_subreads.fasta $output/AAA168-E21
# python $src/coassemble_shuffled.py $options AAA164-B23 $illumina/AAA164B23.ACTTGA/AAA164B23.ACTTGA.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/PB0833_filtered_subreads.fasta $output/AAA164-B23
# python $src/coassemble_shuffled.py $options AAA164-I21 $illumina/AAA164I21.TGACCA/AAA164I21.TGACCA.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/PB0834_filtered_subreads.fasta $output/AAA164-I21
# python $src/coassemble_shuffled.py $options AAA160-C11 $illumina/AAA160C11.GCCAAT/AAA160C11.GCCAAT.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/PB0835_filtered_subreads.fasta $output/AAA160-C11
# python $src/coassemble_shuffled.py $options AAA164-A21 $illumina/AAA164A21.ACTTGA/AAA164A21.ACTTGA.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/PB0836_filtered_subreads.fasta $output/AAA164-A21
#
# pacbio=/mnt/scgc/data/PacBio/PacBio_reads/batch1
# python $src/coassemble_shuffled.py $options AAA164M04 $illumina/AAA164M04.CGATGT/AAA164M04.CGATGT.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/AAA164M04.trimmed.fasta $output/AAA164M04
# python $src/coassemble_shuffled.py $options AAA164P11 $illumina/AAA164P11.GATCAG/AAA164P11.GATCAG.kmernorm.lowcomp_filter.0.010000.fastq.gz $pacbio/AAA164P11.trimmed.fasta $output/AAA164P11

illumina=/mnt/scgc/data/PacBio/Illumina_data
output=/mnt/scgc/data/PacBio/Assemblies/illumina_no_filter
for fastq in $illumina/*.fastq; do
    bn=$(basename $fastq)
    sample=${bn/.fastq}
    echo "processing $sample"
    python $src/assemble_shuffled.py $options $fastq $output/$sample
done
