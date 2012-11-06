#! /usr/bin/env bash
#BSUB -J align
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

BOWTIEINDEX=/vol1/home/brownj/projects/ref/snRNA/snRNA
DATA="$HOME/projects/zhao/data"
FASTA="$DATA/snRNA.fa"
FASTQ="$DATA/s_6_sequence.txt.barcode-TG.fastq"

GAPOPEN=5
GAPEXTEND=3

OPTIONS="-q --rdg $GAPOPEN,$GAPEXTEND -M 10 -R 10 -L 10 -i S,1,1.25 -N 1 -p 8 -x $BOWTIEINDEX -U $FASTQ"

SAMPLE="s_6_sequence.TG.$GAPOPEN.$GAPEXTEND.bowtie2"
UBAM="$SAMPLE.unsorted.bam"

#alignment
bowtie2 $OPTIONS | samtools view -Sb - > $UBAM

#variant calling requires sort
samtools sort $UBAM $SAMPLE

#variant calling
samtools mpileup -ugf $FASTA "$SAMPLE.bam" | bcftools view -vcg - > "$SAMPLE.vcf"

#vcf to bedgraph
python /vol1/home/brownj/projects/hits-clip/src/bcf-identify-cims.py -i "$SAMPLE.vcf" \
    > "$SAMPLE.cims.bedgraph"