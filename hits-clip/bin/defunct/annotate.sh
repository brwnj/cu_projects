#!/usr/bin/env bash
#BSUB -J mapping_summaries[1-47]
#BSUB -e ms.%J.%I.err
#BSUB -o ms.%J.%I.out
#BSUB -q idle
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"

<<DOC
summarize the reads before and after removing duplicates.

need to add snoRNA
DOC

set -o nounset -o pipefail -o errexit -x

samples=(idx0 MP10 MP1 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
           MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG MP42.TCGA MP43.ACTG 
           MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11 
           PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53 
           PK54 helaa helab)
sample=${samples[${LSB_JOBINDEX}]}

intron=$HOME/ref/hg18/refseq.intron.bed.gz
utr3=$HOME/ref/hg18/refseq.utr3.bed.gz
utr5=$HOME/ref/hg18/refseq.utr5.bed.gz
cds=$HOME/ref/hg18/refseq.cds.bed.gz
exon=$HOME/ref/hg18/refseq.exon.bed.gz
whole_gene=$HOME/ref/hg18/refseq.wholegene.bed.gz

mirna_intronic=$HOME/projects/hits-clip/data/common/mirbase/mirbase-18/hsa.hg18.intronic.bed.gz
mirna_intergenic=$HOME/projects/hits-clip/data/common/mirbase/mirbase-18/hsa.hg18.intergenic.bed.gz

datadir=$HOME/projects/hits-clip/results/common/samples/$sample
fullbam=$datadir/$sample.bam
fullbed=$datadir/$sample.bed.gz
smallbam=$datadir/$sample.rmd.bam
smallbed=$datadir/$sample.rmd.bed.gz

if [ ! -f $fullbed ]; then
    bedtools bamtobed -i $fullbam | gzip -c > $fullbed
fi
# read must overlap at least 51%. strand must match. print a once if any overlap.
bedtools intersect -s -f 0.51 -u -a $fullbed -b $intron > $datadir/intron.o.bed
bedtools intersect -s -f 0.51 -u -a $fullbed -b $utr3 > $datadir/utr3.o.bed
bedtools intersect -s -f 0.51 -u -a $fullbed -b $utr5 > $datadir/utr5.o.bed
bedtools intersect -s -f 0.51 -u -a $fullbed -b $cds > $datadir/cds.o.bed
bedtools intersect -s -f 0.51 -u -a $fullbed -b $exon > $datadir/exon.o.bed
bedtools intersect -s -f 0.51 -u -a $fullbed -b $whole_gene > $datadir/whole.o.bed
bedtools intersect -s -f 0.10 -u -a $fullbed -b $mirna_intronic > $datadir/mirna_intronic.o.bed
bedtools intersect -s -f 0.10 -u -a $fullbed -b $mirna_intergenic > $datadir/mirna_intergenic.o.bed

if [ ! -f $smallbed ]; then
    bedtools bamtobed -i $smallbam | gzip -c > $smallbed
fi

bedtools intersect -s -f 0.51 -u -a $smallbed -b $intron > $datadir/intron.d.bed
bedtools intersect -s -f 0.51 -u -a $smallbed -b $utr3 > $datadir/utr3.d.bed
bedtools intersect -s -f 0.51 -u -a $smallbed -b $utr5 > $datadir/utr5.d.bed
bedtools intersect -s -f 0.51 -u -a $smallbed -b $cds > $datadir/cds.d.bed
bedtools intersect -s -f 0.51 -u -a $smallbed -b $exon > $datadir/exon.d.bed
bedtools intersect -s -f 0.51 -u -a $smallbed -b $whole_gene > $datadir/whole.d.bed
bedtools intersect -s -f 0.10 -u -a $smallbed -b $mirna_intronic > $datadir/mirna_intronic.d.bed
bedtools intersect -s -f 0.10 -u -a $smallbed -b $mirna_intergenic > $datadir/mirna_intergenic.d.bed

echo "full counts: intron, utr3, utr5, cds, exon, whole, intronic, intergenic"
wc -l $datadir/intron.o.bed
wc -l $datadir/utr3.o.bed
wc -l $datadir/utr5.o.bed
wc -l $datadir/cds.o.bed
wc -l $datadir/exon.o.bed
wc -l $datadir/whole.o.bed
wc -l $datadir/mirna_intronic.o.bed
wc -l $datadir/mirna_intergenic.o.bed
echo "full bam"
samtools view -cF4 $fullbam
echo "small counts: intron, utr3, utr5, cds, exon, whole, intronic, intergenic"
wc -l $datadir/intron.d.bed
wc -l $datadir/utr3.d.bed
wc -l $datadir/utr5.d.bed
wc -l $datadir/cds.d.bed
wc -l $datadir/exon.d.bed
wc -l $datadir/whole.d.bed
wc -l $datadir/mirna_intronic.d.bed
wc -l $datadir/mirna_intergenic.d.bed
echo "small bam"
samtools view -cF4 $smallbam