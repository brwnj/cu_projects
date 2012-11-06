#! /usr/bin/env bash
#BSUB -J clip[1-2]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

<<DOC

samtools rmdup -s helaa.bam helaa.rmd.bam
bedtools bamtobed -i helaa.bam | gzip -c > helaa.bed.gz
bedtools bamtobed -i helaa.rmd.bam | gzip -c > helaa.rmd.bed.gz
zcat helaa.rmd.bed.gz | awk '$6=="+"' > helaa.rmd.pos.bed
zcat helaa.rmd.bed.gz | awk '$6=="-"' > helaa.rmd.neg.bed
zcat helaa.bed.gz | awk '$6=="+"' > helaa.pos.bed
zcat helaa.bed.gz | awk '$6=="-"' > helaa.neg.bed
bedSort helaa.pos.bed helaa.pos.bed
bedSort helaa.neg.bed helaa.neg.bed 
bedSort helaa.rmd.pos.bed helaa.rmd.pos.bed
bedSort helaa.rmd.neg.bed helaa.rmd.neg.bed 
bedItemOverlapCount -chromSize=/vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes hg18 stdin < helaa.pos.bed > helaa.pos.bedgraph
bedItemOverlapCount -chromSize=/vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes hg18 stdin < helaa.neg.bed > helaa.neg.bedgraph
bedItemOverlapCount -chromSize=/vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes hg18 stdin < helaa.rmd.pos.bed > helaa.rmd.pos.bedgraph
bedItemOverlapCount -chromSize=/vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes hg18 stdin < helaa.rmd.neg.bed > helaa.rmd.neg.bedgraph
bedGraphToBigWig helaa.pos.bedgraph /vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes helaa.pos.bw
bedGraphToBigWig helaa.neg.bedgraph /vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes helaa.neg.bw
bedGraphToBigWig helaa.rmd.pos.bedgraph /vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes helaa.rmd.pos.bw
bedGraphToBigWig helaa.rmd.neg.bedgraph /vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes helaa.rmd.neg.bw

genomedata-open-data combined.genomedata helaa.pos helaa.neg helab.pos helab.neg
zcat ../results/common/samples/helaa/helaa.pos.bedgraph.gz | genomedata-load-data combined.genomedata helaa.pos
zcat ../results/common/samples/helaa/helaa.neg.bedgraph.gz | genomedata-load-data combined.genomedata helaa.neg
zcat ../results/common/samples/helab/helab.pos.bedgraph.gz | genomedata-load-data combined.genomedata helab.pos
zcat ../results/common/samples/helab/helab.neg.bedgraph.gz | genomedata-load-data combined.genomedata helab.neg
genomedata-close-data combined.genomedata
genomedata-open-data combined.rmd.genomedata helaa.pos.rmd helaa.neg.rmd helab.pos.rmd helab.neg.rmd
zcat ../results/common/samples/helaa/helaa.pos.rmd.bedgraph.gz | genomedata-load-data combined.rmd.genomedata helaa.pos.rmd
zcat ../results/common/samples/helaa/helaa.neg.rmd.bedgraph.gz | genomedata-load-data combined.rmd.genomedata helaa.neg.rmd
zcat ../results/common/samples/helab/helab.pos.rmd.bedgraph.gz | genomedata-load-data combined.rmd.genomedata helab.pos.rmd
zcat ../results/common/samples/helab/helab.neg.rmd.bedgraph.gz | genomedata-load-data combined.rmd.genomedata helab.neg.rmd
genomedata-close-data combined.rmd.genomedata
DOC

# Samples for the job array
SAMPLEIDS=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
            MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG MP42.TCGA MP43.ACTG 
            MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11 
            PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53 
            PK54)
SAMPLEIDS=(helaa helab)
SAMPLE_IDX=$(expr $LSB_JOBINDEX - 1)
# Individual sample
SAMPLE=${SAMPLEIDS[$SAMPLE_IDX]}

BCFTOBEDGRAPH=/vol1/home/brownj/projects/hits-clip/src/bcf-identify-cims.py
FASTA=/vol3/home/jhessel/projects/encode/data/hg18/fasta/combined/hg18.fa
CHROMSIZES=/vol3/home/jhessel/projects/encode/data/hg18/hg18.chrom.sizes

# All links created in advance
RAWDATA="/vol1/home/brownj/projects/hits-clip/data/common/$(echo $SAMPLE | cut -d'.' -f1)"
RUNDIR="/vol1/home/brownj/projects/hits-clip/results/common/samples/$SAMPLE"
mkdir $RUNDIR

BAM=$SAMPLE.bam
BED=$SAMPLE.bed
# rmdups
RBAM=$SAMPLE.rmd.bam
RBED=$SAMPLE.rmd.bed

# writes align.sh, stderr, and stdout in the sample's directory
RUNSCRIPT="$RUNDIR/align.sh"
if [ -f $RUNSCRIPT ]; then
    echo "Found: $RUNSCRIPT
          Skipping..."
else
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J $SAMPLE.novoalign" >> $RUNSCRIPT
    echo "#BSUB -n 1" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1]\"" >> $RUNSCRIPT
    echo "#BSUB -e %J.err" >> $RUNSCRIPT
    echo "
NOVOALIGN=$(which novoalign)
NOVOINDEX=/vol1/home/brownj/projects/hits-clip/data/common/novoalign/hg18
SAMTOOLS=$(which samtools)
FASTQ=$RAWDATA/$SAMPLE.fastq.gz
gunzip \$FASTQ
FASTQ=$RAWDATA/$SAMPLE.fastq
OPTIONS=\"-d \$NOVOINDEX -f \$FASTQ -a -o SAM -r A 20 -e 100 -s 2 -l 16\"

\$NOVOALIGN \$OPTIONS | \$SAMTOOLS view -Sb - > $SAMPLE.unsorted.bam
\$SAMTOOLS sort $SAMPLE.unsorted.bam $SAMPLE
rm -f $SAMPLE.unsorted.bam
" >> $RUNSCRIPT

    # Submit alignment
    bsub -cwd $RUNDIR < $RUNSCRIPT
fi

RUNSCRIPT="$RUNDIR/collapse.sh"
if [ -f $RUNSCRIPT ]; then
    echo "Found: $RUNSCRIPT
          Skipping..."
else
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J $SAMPLE.collapse" >> $RUNSCRIPT
    echo "#BSUB -o %J.out" >> $RUNSCRIPT
    echo "#BSUB -e %J.err" >> $RUNSCRIPT
    echo "#BSUB -n 1" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1]\"" >> $RUNSCRIPT
    echo "
samtools rmdup -s $BAM $RBAM
bedtools bamtobed -i $RBAM | gzip -c > $RBED
" >> $RUNSCRIPT
    WAIT=$(bjobs -w | grep $SAMPLE.align | awk '{print $1}' | sort | uniq)
    if [ x$WAIT == "x" ]; then
        bsub -cwd $RUNDIR < $RUNSCRIPT 
    else
        bsub -w $WAIT -cwd $RUNDIR < $RUNSCRIPT
    fi
fi

# Start data summary
RUNSCRIPT="$RUNDIR/summary.sh"
if [ -f $RUNSCRIPT ]; then
    echo "Found: $RUNSCRIPT
          Skipping..."
else
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J $SAMPLE.summary" >> $RUNSCRIPT
    echo "#BSUB -o %J.out" >> $RUNSCRIPT
    echo "#BSUB -e %J.err" >> $RUNSCRIPT
    echo "#BSUB -n 1" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1]\"" >> $RUNSCRIPT
    echo "
COMBINEDBED=$RBED
POSBED=$SAMPLE.pos.rmd.bed
NEGBED=$SAMPLE.neg.rmd.bed
POSBG=$SAMPLE.pos.rmd.bedgraph
NEGBG=$SAMPLE.neg.rmd.bedgraph
POSBW=$SAMPLE.pos.rmd.bw
NEGBW=$SAMPLE.neg.rmd.bw

# Create a bigwig for pos and neg strand
zcat \$COMBINEDBED | awk '\$6==\"+\"' > \$POSBED
zcat \$COMBINEDBED | awk '\$6==\"-\"' > \$NEGBED
bedSort \$POSBED \$POSBED
bedSort \$NEGBED \$NEGBED
bedItemOverlapCount -chromSize=$CHROMSIZES hg18 stdin < \$POSBED > \$POSBG
bedItemOverlapCount -chromSize=$CHROMSIZES hg18 stdin < \$NEGBED > \$NEGBG
bedGraphToBigWig \$POSBG $CHROMSIZES \$POSBW
bedGraphToBigWig \$NEGBG $CHROMSIZES \$NEGBW

gzip *.bed
gzip *.bedgraph

# Create bedgraph for alignment
bedtools genomecov -bg -ibam $RBAM -g $CHROMSIZES | gzip -c > $SAMPLE.rmd.bedgraph.gz
" >> $RUNSCRIPT

    # Have this job submitted, but wait on the alignment
    WAIT=$(bjobs -w | grep $SAMPLE.collapse | awk '{print $1}' | sort | uniq)
    if [ x$WAIT == "x" ]; then
        bsub -q normal -cwd $RUNDIR < $RUNSCRIPT 
    else
        bsub -q normal -w $WAIT -cwd $RUNDIR < $RUNSCRIPT
    fi
fi

# Start CIMS calling
# RUNSCRIPT="$RUNDIR/cims.sh"
# if [ -f $RUNSCRIPT ]; then
#     echo "Found: $RUNSCRIPT
#           Skipping..."
# else
#     echo "#! /usr/bin/env bash" > $RUNSCRIPT
#     echo "#BSUB -J $SAMPLE.CIMS" >> $RUNSCRIPT
#     echo "#BSUB -n 1" >> $RUNSCRIPT
#     echo "#BSUB -R \"span[hosts=1]\"" >> $RUNSCRIPT
#     echo "#BSUB -e %J.err" >> $RUNSCRIPT
#     echo "#BSUB -o %J.out" >> $RUNSCRIPT
#     echo "
# samtools mpileup -6 -ugf $FASTA $SAMPLE.bam | bcftools view -bvcg - > $SAMPLE.bcf
# bcftools view $SAMPLE.bcf | python $BCFTOBEDGRAPH -i - > $SAMPLE.novoalign.cims.bedgraph
# bedGraphToBigWig $SAMPLE.novoalign.cims.bedgraph $CHROMSIZES $SAMPLE.novoalign.cims.bw
# # Transfer bigwig to sandbox
# # scp *.bw amc-sandbox:/data/home/brownj/public_html/hits-clip/
# " >> $RUNSCRIPT
#     # Again, submit but wait on alignment
#     if [ x$WAIT == "x" ]; then
#         bsub -q normal -cwd $RUNDIR < $RUNSCRIPT 
#     else
#         bsub -q normal -w $WAIT -cwd $RUNDIR < $RUNSCRIPT
#     fi
# fi

# Write no matter what to create a new session
# Write UCSC tracks for alignment and CIMS
# TRACKS=tracks.txt
# echo "track type=bigWig name='$SAMPLE.rmd.pos.align' description='$SAMPLE Novoalign alignment pos ($(date +%Y%m%d))' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/hits-clip/$SAMPLE.pos.rmd.bw maxHeightPixels=15:50:35 color=252,141,98 visibility=full" >> $TRACKS
# echo "track type=bigWig name='$SAMPLE.rmd.neg.align' description='$SAMPLE Novoalign alignment neg ($(date +%Y%m%d))' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/hits-clip/$SAMPLE.neg.rmd.bw maxHeightPixels=15:50:35 color=252,141,98 visibility=full" >> $TRACKS
# echo "track type=bigWig name='$SAMPLE.novo.CIMS' description='$SAMPLE Novoalign CIMS ($(date +%Y%m%d))' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/hits-clip/$SAMPLE.novoalign.cims.bw maxHeightPixels=15:50:35 color=252,141,98 visibility=full" >> $TRACKS