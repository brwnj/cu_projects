#! /usr/bin/env bash
#BSUB -J cleanup_peaks
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q short

<<DOC
Cleanup the shell scripts and concat chromosome files.
DOC

set -o nounset -o pipefail -o errexit -x

# SAMPLES=(idx0 MCF7 MCF7estr MDA231 BT474 BT474estr BT474herc HS27A HS5 hMSC BMEC HUVEC)
# SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}
SAMPLE=HELA

CHROMSIZES=$HOME/ref/hg18/chrom_sizes.txt

DATA=$HOME/projects/hits-clip/results/common
RUNDIR=$DATA/$SAMPLE

RUN=$RUNDIR/$SAMPLE.cleanup.sh

if [ ! -f $RUNDIR/$SAMPLE.pos.peaks.trimmed.bed.gz -a ! -f $RUNDIR/$SAMPLE.neg.peaks.trimmed.bed.gz ]; then
    echo "#! /usr/bin/env bash" > $RUN
    echo "#BSUB -J $SAMPLE.cleanup" >> $RUN
    echo "#BSUB -o %J.out" >> $RUN
    echo "#BSUB -e %J.err" >> $RUN
    echo "#BSUB -q short" >> $RUN
    echo "
    STRANDS=\"pos neg\"
    CHROMOSOMES=\"chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11
                 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21
                 chr22 chrX chrY chrM\"

    for STRAND in \$STRANDS; do
        COMBINEDBED=$SAMPLE.\$STRAND.peaks.trimmed.combined.bed

        for CHROM in \$CHROMOSOMES; do
            BED=$SAMPLE.\$CHROM.\$STRAND.peaks.trimmed.bed.gz
            zcat \$BED >> \$COMBINEDBED
            rm -f \$BED
        done
        bedSort \$COMBINEDBED \$COMBINEDBED

        FIXED=$SAMPLE.\$STRAND.peaks.trimmed.fixed.bed
        CLIPPED=$SAMPLE.\$STRAND.peaks.trimmed.bed
        TRACK=$SAMPLE.\$STRAND.peaks.trimmed.track.bed
        bioawk -c header '\$2<\$3' \$COMBINEDBED > \$FIXED
        bedClip -verbose=2 \$FIXED $CHROMSIZES \$CLIPPED
        echo \"track type=bed name=$SAMPLE.\$STRAND.60bp visibility=pack\" > \$TRACK
        cat \$CLIPPED >> \$TRACK
        gzip -f \$TRACK
        gzip -f \$CLIPPED
        rm -f \$COMBINEDBED \$FIXED
    done

    rm -f *trimpeaks.sh

    " >> $RUN

    bsub -cwd $RUNDIR < $RUN
fi