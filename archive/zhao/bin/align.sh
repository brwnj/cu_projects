#! /usr/bin/env bash
#BSUB -J align
#BSUB -e %J.err
#BSUB -o %J.out
#BSUB -q normal

GAPOPENINGPENALTY=10
GAPEXTENDPENALTY=1

NOVOALIGN="$(which novoalign)"
DATA="$HOME/projects/zhao/data"
NOVOINDEX="$DATA/snRNA"

FASTA="$DATA/snRNA.fa"

FASTQ="$DATA/s_6_sequence.txt.barcode-TG.fastq"
OPTIONS="-d $NOVOINDEX -f $FASTQ -o SAM -r N -c 12 -g $GAPOPENINGPENALTY -x $GAPEXTENDPENALTY"

SAMPLE="s_6_sequence.TG.$GAPOPENINGPENALTY.$GAPEXTENDPENALTY"
UBAM="$SAMPLE.unsorted.bam"

#alignment
$NOVOALIGN $OPTIONS | samtools view -Sb - > $UBAM

#variant calling requires sort
samtools sort $UBAM $SAMPLE

#variant calling
samtools mpileup -ugf $FASTA "$SAMPLE.bam" | bcftools view -vcg - > "$SAMPLE.vcf"

#vcf to bedgraph
python /vol1/home/brownj/projects/hits-clip/src/bcf-identify-cims.py -i "$SAMPLE.vcf" \
    > "$SAMPLE.cims.bedgraph"