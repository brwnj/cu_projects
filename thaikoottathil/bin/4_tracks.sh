#! /usr/bin/env bash
#BSUB -J ucsc.tracks[1-13]
#BSUB -e tracks.%J.%I.err
#BSUB -o tracks.%J.%I.out
#BSUB -q normal

<<DOC
generate bigwigs of the alignment and write track lines to a stdout.
DOC

PI=thaikoottathil

set -o nounset -o pipefail -o errexit -x

SAMPLES=(idx
10_AAGCTA_L006_R1_001.fastq
11_GTAGCC_L006_R1_001.fastq
12_TACAAG_L006_R1_001.fastq
14_GGAACT_L006_R1_001.fastq
15_TGACAT_L006_R1_001.fastq
16_GGACGG_L006_R1_001.fastq
18_TTTCAC_L006_R1_001.fastq
19_GGCCAC_L006_R1_001.fastq
22_TCAAGT_L006_R1_001.fastq
23_ACATCG_L006_R1_001.fastq
24_ATTGGC_L006_R1_001.fastq
3_GCCTAA_L006_R1_001.fastq
7_GATCTG_L006_R1_001.fastq
)

SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}
NAME=$(basename $SAMPLE .fastq)
PROJECT=$HOME/projects/$PI

# UCSC Browser stuff
COMBINEDBED=$NAME.combined.bed.gz
POSBED=$NAME.pos.bed
NEGBED=$NAME.neg.bed
POSBG=$NAME.pos.bedgraph
NEGBG=$NAME.neg.bedgraph
POSBW=$NAME.pos.bw
NEGBW=$NAME.neg.bw
CHROMSIZES=/vol1/home/brownj/projects/ref/hg19/chrom_sizes.txt

if [ ! -f $POSBW ]; then
    # Create a bigwig for pos and neg strand
    bedtools bamtobed -i $PROJECT/results/common/$NAME/$NAME.bam | awk 'length($1)<6 && $1!="chrM"' | gzip -c > $COMBINEDBED
    zcat $COMBINEDBED | awk '$6=="+"' > $POSBED
    zcat $COMBINEDBED | awk '$6=="-"' > $NEGBED
    bedSort $POSBED $POSBED
    bedSort $NEGBED $NEGBED
    bedItemOverlapCount -chromSize=$CHROMSIZES hg19 stdin < $POSBED > $POSBG
    bedItemOverlapCount -chromSize=$CHROMSIZES hg19 stdin < $NEGBED > $NEGBG
    bedGraphToBigWig $POSBG $CHROMSIZES $POSBW
    bedGraphToBigWig $NEGBG $CHROMSIZES $NEGBW

    # CLean up
    gzip $POSBED $NEGBED
    gzip $POSBG $NEGBG
fi

TRACK=$NAME.tracks.txt

echo "track type=bed name=${NAME}_peaks visibility=pack" > $NAME.peaks.track.bed
zcat $PROJECT/results/common/$NAME/${NAME}_peaks.bed.gz >> $NAME.peaks.track.bed
gzip $NAME.peaks.track.bed

echo "track type=bed name=${NAME}_subpeaks visibility=pack" > $NAME.subpeaks.track.bed
zcat $PROJECT/results/common/$NAME/${NAME}_peaks.subpeaks.bed.gz | awk '$1!="Chromosome"' >> $NAME.subpeaks.track.bed
gzip $NAME.subpeaks.track.bed

echo "track type=bigWig name='$NAME pos' description='$NAME pos' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/$PI/$POSBW maxHeightPixels=15:50:35 color=77,175,74 visibility=full" > $TRACK
echo "track type=bigWig name='$NAME neg' description='$NAME neg' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/$PI/$NEGBW maxHeightPixels=15:50:35 color=77,175,74 visibility=full" >> $TRACK

echo "http://amc-sandbox.ucdenver.edu/~brownj/$PI/$NAME.peaks.track.bed.gz" >> $TRACK
echo "http://amc-sandbox.ucdenver.edu/~brownj/$PI/$NAME.subpeaks.track.bed.gz" >> $TRACK

echo "http://amc-sandbox.ucdenver.edu/~brownj/$PI/${NAME}_treat_afterfiting_all.wig.gz" >> $TRACK
echo "http://amc-sandbox.ucdenver.edu/~brownj/$PI/${NAME}_control_afterfiting_all.wig.gz" >> $TRACK