#! /usr/bin/env bash
#BSUB -J trim_peaks
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal

<<DOC
of the combined peaks, find the midpoints of each peak, and expand to the size
of the ago footprint.
DOC
# 
# SAMPLES=(MCF7 MCF7estr MDA231 BT474 BT474estr BT474herc HS27A HS5 hMSC BMEC HUVEC)
# REPLICATES=("PK11 PK31 PK51" "PK12 PK32" "PK24 PK42 PK54" "PK21 PK41 PK52" \
#             "PK22 PK53" "PK23" "MP1 MP21 MP35" "MP2 MP20 MP34"\
#             "MP36 MP43.ACTG MP43.TCGA MP44.ACTG MP44.TCGA"\
#             "MP42.ACTG MP45.ACTG MP45.TCGA" "MP24 MP38")
# IDX=$(expr $LSB_JOBINDEX - 1)
# SAMPLE=${SAMPLES[$IDX]}
# REPLICATE=${REPLICATES[$IDX]}
SAMPLE=HUVEC
REPLICATE="MP24 MP38"

DATA=/vol1/home/brownj/projects/hits-clip/results/common
RUNDIR=$DATA/$SAMPLE

FOOTPRINT=60
SRC=$HOME/projects/hits-clip/bin
GENOMEDATA=$PROJECT/data/combined.rmd.genomedata
#GENOMEDATA=/vol1/home/brownj/projects/hits-clip/results/20120917/gd_20120917.rmd

STRANDS="pos neg"
CHROMOSOMES="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 \
             chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 \
             chr22 chrX chrY chrM"

for STRAND in $STRANDS; do        
    if [ ! -f $RUNDIR/$SAMPLE.$STRAND.peaks.trimmed.bed.gz ]; then  
        for CHROM in $CHROMOSOMES; do
            BED=$SAMPLE.$STRAND.peaks.bed.gz
            TRIMMEDPEAKS=$SAMPLE.$CHROM.$STRAND.peaks.trimmed.bed.gz
            TRACKS="$(for R in $REPLICATE; do echo -n "$R.$STRAND.rmd "; done)"
            OPTIONS="-w $FOOTPRINT -c $CHROM $BED $GENOMEDATA $TRACKS"
            RUN=$RUNDIR/$STRAND.$CHROM.trimpeaks.sh
        
            echo "#! /usr/bin/env bash" > $RUN
            echo "#BSUB -J $SAMPLE.$STRAND.$CHROM.trimming" >> $RUN
            echo "#BSUB -o %J.out" >> $RUN
            echo "#BSUB -e %J.err" >> $RUN
            echo "" >> $RUN
            echo "python $SRC/trim_peaks.py $OPTIONS | gzip -c > $TRIMMEDPEAKS" >> $RUN
        
            bsub -cwd $RUNDIR < $RUN
        done
    fi
done