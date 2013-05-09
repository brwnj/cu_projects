#! /usr/bin/env bash
#BSUB -J peaktracks
#BSUB -o %J.out
#BSUB -e %J.err
#BSUB -q short

<<DOC
Creates the necessary files for the individual sample tracks to
compare alongside the combined replicate files.
DOC

# SAMPLEIDS=(MP1 MP10 MP11 MP2 MP20 MP21 MP22 MP23 MP24 MP30 MP31 MP34 MP35 MP36
#            MP38 MP39.ACTG MP39.TCGA MP40 MP41 MP42.ACTG MP42.TCGA MP43.ACTG 
#            MP43.TCGA MP44.ACTG MP44.TCGA MP45.ACTG MP45.TCGA MP7 MP9 PK11 
#            PK12 PK21 PK22 PK23 PK24 PK31 PK32 PK33 PK41 PK42 PK51 PK52 PK53 
#            PK54)
# IDX=$(expr $LSB_JOBINDEX - 1)
# SAMPLE=${SAMPLEIDS[$IDX]}

CHROMSIZES=$HOME/ref/hg18/hg18.sizes
CUTOFF=.001

BED="$SAMPLE.peaks.rmd.qv$CUTOFF.bed.gz"
FIXEDBED="$SAMPLE.peaks.rmd.f.bed"
CLIPPEDBED="$SAMPLE.peaks.rmd.f.c.bed"
TRACKBED="$SAMPLE.peaks.track$CUTOFF.bed"

RUNDIR=$HOME/projects/hits-clip/results/common/samples/$SAMPLE
RUN=$RUNDIR/peaktracks.sh
echo "#!/usr/bin/env bash" > $RUN
echo "#BSUB -J $SAMPLE.bedtrack" >> $RUN
echo "#BSUB -o %J.out" >> $RUN
echo "#BSUB -e %J.err" >> $RUN
echo "
awk -c header '\$2<\$3{print \$1,\$2,\$3,\$7,\$7*1000,\$6}' $BED > $FIXEDBED
bedClip -verbose=2 $FIXEDBED $CHROMSIZES $CLIPPEDBED
echo 'track type=bed name=$SAMPLE.peaks visibility=pack useScore=1' > $TRACKBED
cat $CLIPPEDBED >> $TRACKBED
gzip -f $TRACKBED
rm -f $FIXEDBED $BED
gzip $CLIPPEDBED
mv $CLIPPEDBED.gz $BED
" >> $RUN
bsub -cwd $RUNDIR < $RUN
