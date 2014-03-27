#!/usr/bin/env bash
#BSUB -J makehub
#BSUB -e makehub.%J.err
#BSUB -o makehub.%J.out
#BSUB -q normal
#BSUB -n 1
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh


# genomes.txt
if [[ ! -f $HUB/genomes.txt ]]; then
    genomes=$HUB/genomes.txt
    echo "genome $GENOME" > $genomes
    echo "trackDb $GENOME/trackDb.txt" >> $genomes
fi

# hub.txt
if [[ ! -f $HUB/hub.txt ]]; then
    hub=$HUB/hub.txt
    echo "hub williamschip" > $hub
    echo "shortLabel Williams ChIP" >> $hub
    echo "longLabel Williams ChIP" >> $hub
    echo "genomesFile genomes.txt" >> $hub
    echo "email brwnjm@gmail.com" >> $hub
fi

# data dir
if [[ ! -d $HUB/$GENOME ]]; then
    mkdir -p $HUB/$GENOME
fi

trackdb=$HUB/$GENOME/trackDb.txt

# output for coverage
cat <<coverage_track >$trackdb
track williams_coverage
compositeTrack on
shortLabel Coverage
longLabel Coverage
maxHeightPixels 50:20:15
subGroup1 cline CellLine s3B5_hela=s3B5_hela s3B5_mitotic_hela=s3B5_mitotic_hela polII_hela=polII_hela polII_mitotic_hela=polII_mitotic_hela SC184_hela=SC184_hela SC184_mitotic_hela=SC184_mitotic_hela hela=hela mitotic_hela=mitotic_hela
type bigWig
configurable on
autoScale on

coverage_track

for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    cp $RESULTS/$sample/*.bw $HUB/$GENOME
    posbw=${sample}_pos.bw
    negbw=${sample}_neg.bw
    cat <<coverage_track >>$trackdb
    track ${posbw/.bw}
    bigDataUrl $posbw
    shortLabel $sample coverage POS
    longLabel $sample coverage positive (+) strand
    subGroups cline=${UCSC_GROUP[$sample]}
    type bigWig
    parent williams_coverage
    color ${COLORS[$sample]}

    track ${negbw/.bw}
    bigDataUrl $negbw
    shortLabel $sample coverage NEG
    longLabel $sample coverage negative (-) strand
    subGroups cline=${UCSC_GROUP[$sample]}
    type bigWig
    parent williams_coverage
    color ${COLORS[$sample]}

coverage_track
done

# output for peaks
cat <<peak_track >>$trackdb
track williams_peaks
compositeTrack on
configurable on
shortLabel Peaks
longLabel Peaks
subGroup1 cline CellLine s3B5_hela=s3B5_hela s3B5_mitotic_hela=s3B5_mitotic_hela polII_hela=polII_hela polII_mitotic_hela=polII_mitotic_hela SC184_hela=SC184_hela SC184_mitotic_hela=SC184_mitotic_hela hela=hela mitotic_hela=mitotic_hela
type bed 6 .

peak_track

for bigbed in $RESULTS/*/*peaks.bb; do
    cp $bigbed $HUB/$GENOME
    bb=$(basename $bigbed)
    sample=${bb/_peaks.bb}
    cat <<peak_track >>$trackdb
    track ${bb/.bb}
    parent williams_peaks
    shortLabel ${bb/_peaks.bb} peaks
    longLabel ${bb/_peaks.bb} peaks
    subGroups cline=${UCSC_GROUP[$sample]}
    bigDataUrl $bb
    color ${COLORS[$sample]}
    thickDrawItem on
    type bigBed 6 .

peak_track
done
