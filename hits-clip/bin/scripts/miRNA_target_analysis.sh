#! /usr/bin/env bash
#BSUB -J "mirna_targets[1-29]%8"
#BSUB -o mirna_targets.%J.%I.out
#BSUB -e mirna_targets.%J.%I.err
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 8
#BSUB -P hits-clip

set -o nounset -o errexit -o pipefail -x

<<DOC
align to mirbase,
    convert alignment to fastq,
    align to genome -r None,
    mask fastq, then
    realign to peak regions.
later, run mirza on the output.
DOC

sampleids=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP34 MP35 MP36 MP7
            MP9 PK11 PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51
            PK52 PK53 PK54) # PK61 PK62 PK63)
id=${sampleids[$(($LSB_JOBINDEX - 1))]}
bin=$HOME/devel/cliptools/cliptools

FASTQ=$HOME/projects/hits-clip/data/common/$id/$id.fastq.gz
MIRBASE_IDX=/vol3/home/jhessel/projects/hits-clip/data/20130121/mature.hsa.mirbase19.ndx
MIRBASE_ALIGN=$id.mirbase19.bam

MIRBASE_FASTQ=$id.mirbase19.fastq.gz

GENOME_IDX=$HOME/projects/hits-clip/data/common/novoalign/hg18
GENOME_ALIGN=$id.genome.bam

MASK_FASTQ=$id.masked.fastq.gz

PEAKS_ALIGN=$id.peaks.bam
PEAKS_IDX=/vol3/home/jhessel/projects/hits-clip/data/20130130/$id.peaks.nidx

# align to mirbase
if [[ ! -f $MIRBASE_ALIGN ]]; then
    novoalign -d $MIRBASE_IDX -f $FASTQ \
            -a -o SAM -r A -e 100 -s 2 -l 14 -c 8 \
        | samtools view -ShuF4 - \
        | samtools sort -o - $id.temp -m 9500000000 \
        > $MIRBASE_ALIGN
fi
if [[ ! -f $MIRBASE_ALIGN.bai ]]; then
    samtools index $MIRBASE_ALIGN
fi

# convert to fastq and align to genome (-r None)
if [[ ! -f $MIRBASE_FASTQ ]]; then
    java -jar ~/opt/picard-tools-1.79/SamToFastq.jar \
        I=$MIRBASE_ALIGN F=${MIRBASE_FASTQ%.gz} RC=true
    gzip ${MIRBASE_FASTQ%.gz}
fi

# align to whole genome, discarding multiple aligners
if [[ ! -f $GENOME_ALIGN ]]; then
    novoalign -d $GENOME_IDX -f $MIRBASE_FASTQ \
            -a -o SAM -r A -e 100 -s 2 -l 14 -c 8 \
        | samtools view -ShuF4 - \
        | samtools sort -o - $id.temp -m 9500000000 \
        > $GENOME_ALIGN
fi
if [[ ! -f $GENOME_ALIGN.bai ]]; then
    samtools index $GENOME_ALIGN
fi

# mask fastq file
if [[ ! -f $MASK_FASTQ ]]; then
    python $bin/mask_fastq.py \
            --debug $GENOME_ALIGN $FASTQ \
        | gzip -c \
        > $MASK_FASTQ
fi

# realign to hg18 knownGenes
if [[ ! -f $PEAKS_ALIGN ]]; then
    novoalign -d $PEAKS_IDX -f $MASK_FASTQ \
            -a -o SAM -r None -s 2 -l 14 -c 8 \
        | samtools view -ShuF4 - \
        | samtools sort -o - $id.temp -m 9500000000 \
        > $PEAKS_ALIGN
fi
if [[ ! -f $PEAKS_ALIGN.bai ]]; then
    samtools index $PEAKS_ALIGN
fi