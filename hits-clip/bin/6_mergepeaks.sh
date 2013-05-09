#!/usr/bin/env bash
#BSUB -J mergepeaks
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q short

<<DOC
Merges the replicate groups without regard to strand and strand specific. It
also filters the bed where start is greater than stop and clips coordinates
down to chromosome sizes.
DOC

# SAMPLES=(idx0 BT474herc MCF7 MCF7estr MDA231 BT474 BT474estr HS27A HS5 hMSC BMEC HUVEC)
# REPLICATES=(idx0 "PK23" "PK11 PK31 PK51" "PK12 PK32" "PK24 PK42 PK54" \
#             "PK21 PK41 PK52" "PK22 PK53" "MP1 MP21 MP35" "MP2 MP20 MP34"\
#             "MP36 MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA"\
#             "MP42.ACTG MP45.ACTG MP45.TCGA" "MP24 MP38") 
# SAMPLE=${SAMPLES[${LSB_JOBINDEX}]}
# REPLICATE=${REPLICATES[${LSB_JOBINDEX}]}
SAMPLE=something
REPLICATE="PK61 PK62 PK63"
SRC=$HOME/devel/peaktools/peaktools/shit...

CHROMSIZES=$HOME/ref/hg18/hg18.sizes
CASE_DATA=$HOME/projects/hits-clip/results/common
REPLICATE_DATA=$CASE_DATA/samples
EXT=peaks.rmd.qv.001.bed.gz

# combines replicate without respect to strand, positive, and negative.
function combine_replicates()
{
    # sample_name replicates strand
    REPLICATES="$1"
    STRAND=${2:-"both"}
    
    # neither "fixed" nor clipped
    BADBED="$SAMPLE.$STRAND.peaks.bad.bed.gz"
    # exclude cases where start < stop
    FIXED="$SAMPLE.$STRAND.peaks.fixed.bed"
    # remove stuff outside of chromosomal coords
    CLIPPED="$SAMPLE.$STRAND.peaks.bed"
    # UCSC track
    TRACK="$SAMPLE.$STRAND.peaks.track.bed"
    
    if [ ! -f "$REPLICATE_DATA/$SAMPLE/$CLIPPED.gz" ]; then
        if [ x$SAMPLE == "xBT474herc" ]; then
            zcat $REPLICATES | gzip -c > $BADBED
        else
            python $SRC/combine_replicates.py -v $REPLICATES | gzip -c > $BADBED
        fi

        bioawk -c header '$2<$3' $BADBED > $FIXED
        bedClip -verbose=2 $FIXED $CHROMSIZES $CLIPPED
        echo "track type=bed name=$SAMPLE.$STRAND visibility=pack" > $TRACK
        cat $CLIPPED >> $TRACK
        gzip -f $TRACK
        gzip -f $CLIPPED
        rm -f $BADBED $FIXED
    fi
}

# write the full path of each sample replicate and creates pos and neg from 
# qvalue file
function run_replicates()
{
    ARGS="$@"
    POSREPS=""
    NEGREPS=""
    REPS=""
    for R in $ARGS; do
        #before combining, separate into pos and neg
        #ie: PK11.peaks.rmd.qv.001.bed.gz
        PEAK="$REPLICATE_DATA/$R/$R.$EXT"
        POSBED="$REPLICATE_DATA/$R/$R.pos.$EXT"
        NEGBED="$REPLICATE_DATA/$R/$R.neg.$EXT"
        bioawk -c header '$6=="+"' $PEAK | gzip -c > $POSBED
        bioawk -c header '$6=="-"' $PEAK | gzip -c > $NEGBED
        POSREPS="$POSREPS $REPLICATE_DATA/$R/$R.pos.$EXT"
        NEGREPS="$NEGREPS $REPLICATE_DATA/$R/$R.neg.$EXT"
        REPS="$REPS $REPLICATE_DATA/$R/$R.$EXT"
    done
    # strand pos neg and both
    POSRUN=$(combine_replicates "$POSREPS" pos)
    NEGRUN=$(combine_replicates "$NEGREPS" neg)
    BOTHRUN=$(combine_replicates "$REPS")
}

if [ ! -f $CASE_DATA/$SAMPLE/$SAMPLE.pos.peaks.trimmed.bed.gz -a ! -f $CASE_DATA/$SAMPLE/$SAMPLE.neg.peaks.trimmed.bed.gz ]; then
    CONSENSUS=$(run_replicates $REPLICATE)
    # copy stuff over to the common dir
    if [ ! -d $CASE_DATA/$SAMPLE/ ]; then
        mkdir $CASE_DATA/$SAMPLE
    fi
    cp $SAMPLE.* $CASE_DATA/$SAMPLE/
fi
