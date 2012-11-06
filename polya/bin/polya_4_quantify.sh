#! /usr/bin/env bash
#BSUB -J genomecov[1-5]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

<<DOC
use bedtools genome coverage to generate bedgraph of 5' end coverage.
DOC

TODAY=$(date "+%Y%m%d")

SAMPLES=(idx0 MP51 MP52 MP53 PK61 PK62)
SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}
COLORS=(255,255,51 255,127,0 152,78,163 77,175,74 55,126,184 228,26,28)
COLOR=${COLORS[${LSB_JOBINDEX}]}

CHROM=/vol1/home/brownj/projects/ref/hg18/chrom_sizes.txt
BIN=$HOME/projects/hits-clip/bin

MOD=umi.unique

bedtools genomecov -bg -5 -ibam $SAMPLE.$MOD.bam -g $CHROM > $SAMPLE.$MOD.bedgraph
bedGraphToBigWig $SAMPLE.$MOD.bedgraph $CHROM $SAMPLE.$MOD.bw
gzip $SAMPLE.$MOD.bedgraph

echo "track type=bigWig name='$SAMPLE.polya' description='$SAMPLE.polya' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/hits-clip/$TODAY/$SAMPLE.$MOD.bw maxHeightPixels=15:50:35 color=$COLOR visibility=full"