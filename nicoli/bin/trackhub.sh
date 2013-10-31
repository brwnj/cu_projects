#!/usr/bin/env bash
#BSUB -J trackhub
#BSUB -e trackhub.%J.err
#BSUB -o trackhub.%J.out
#BSUB -q normal
#BSUB -n 1
#BSUB -P nicoli

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

genome=hg18

if [[ ! -d $HUB ]]; then
    mkdir $RESULTS/hub
fi
# genomes.txt
if [[ ! -f $HUB/genomes.txt ]]; then
    genomes=$HUB/genomes.txt
    echo "genome $genome" > $genomes
    echo "trackDb $genome/trackDb.txt" >> $genomes
fi
# hub.txt
if [[ ! -f $HUB/hub.txt ]]; then
    hub=$HUB/hub.txt
    echo "hub HITS-CLIP" > $hub
    echo "shortLabel HITS-CLIP" >> $hub
    echo "longLabel HITS-CLIP" >> $hub
    echo "genomesFile genomes.txt" >> $hub
    echo "email brwnjm@gmail.com" >> $hub
fi
# data dir
if [[ ! -d $HUB/$genome ]]; then
    mkdir -p $HUB/$genome
fi

trackdb=$HUB/$genome/trackDb.txt

# output for coverage
cat <<coverage_track >$trackdb
track nicoli_coverage
compositeTrack on
shortLabel Nicoli coverage
longLabel Nicoli coverage
maxHeightPixels 50:20:15
type bigWig
configurable on
autoScale on

coverage_track
for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    cp $RESULTS/$sample/$sample.*.bw $HUB/$genome
    posbw=$sample.rmd_pos.bw
    negbw=$sample.rmd_neg.bw
    color=$(python -c 'import colorbrewer,random;print ",".join(map(str, colorbrewer.Paired[10][random.randint(0,10)]))')
    cat <<coverage_track >>$trackdb
    track ${posbw/.bw}
    bigDataUrl $posbw
    shortLabel $sample coverage POS
    longLabel $sample coverage positive (+) strand
    type bigWig
    parent nicoli_coverage
    color $color
    
    track ${negbw/.bw}
    bigDataUrl $negbw
    shortLabel $sample coverage NEG
    longLabel $sample coverage negative (-) strand
    type bigWig
    parent nicoli_coverage
    color $color

coverage_track
done

# output for peaks
cat <<peak_track >>$trackdb
track nicoli_peaks
compositeTrack on
configurable on
shortLabel Nicoli peaks
longLabel Nicoli peaks (q < 0.001)
type bed 6 .

peak_track
for peaks in $RESULTS/*/*peaks.qv.passed_filter.bed.gz; do
    bb=${peaks/.bed.gz/.bb}
    if [[ ! -f $bb ]]; then
        bed2bb.py --type bed6+1 $SIZES $peaks
    fi
    # check again to make sure previous script didn't fail
    if [[ ! -f $bb ]]; then
        echo "conversion to bb failed on file: $peaks"
        exit 1
    fi
    cp $bb $HUB/$genome
    bb=$(basename $bb)
    cat <<peak_track >>$trackdb
    track ${bb/.bb}
    parent nicoli_peaks
    shortLabel ${bb/.*} peaks
    longLabel ${bb/.*} peaks
    bigDataUrl $bb
    color 31,120,180
    thickDrawItem on
    type bigBed 6 .

peak_track
done
