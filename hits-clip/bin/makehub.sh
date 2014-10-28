#!/usr/bin/env bash
#BSUB -J makehub
#BSUB -e makehub.%J.err
#BSUB -o makehub.%J.out
#BSUB -q normal
#BSUB -n 1
#BSUB -P hits-clip

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/hits-clip/bin/config.sh

# genomes.txt
if [[ ! -f $HUB/genomes.txt ]]; then
    genomes=$HUB/genomes.txt
    echo "genome $GENOME" > $genomes
    echo "trackDb $GENOME/trackDb.txt" >> $genomes
fi

# hub.txt
if [[ ! -f $HUB/hub.txt ]]; then
    hub=$HUB/hub.txt
    echo "hub hitsclip" > $hub
    echo "shortLabel HITS-CLIP" >> $hub
    echo "longLabel HITS-CLIP" >> $hub
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
track hitsclip_coverage
compositeTrack on
shortLabel Coverage
longLabel Coverage
maxHeightPixels 50:20:15
type bigWig
configurable on
autoScale on

coverage_track

for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    posbw=${sample}_pos.bw
    negbw=${sample}_neg.bw
    posuniqbw=$sample.rmd_pos.bw
    neguniqbw=$sample.rmd_neg.bw

    cp $RESULTS/$sample/*.bw $HUB/$GENOME

    color=$(python -c 'import colorbrewer,random;print ",".join(map(str, colorbrewer.Paired[10][random.randint(0,9)]))')
    cat <<coverage_track >>$trackdb
    track ${posbw/.bw}
    bigDataUrl $posbw
    shortLabel $sample coverage POS
    longLabel $sample coverage positive (+) strand
    type bigWig
    parent hitsclip_coverage
    color $color

    track ${negbw/.bw}
    bigDataUrl $negbw
    shortLabel $sample coverage NEG
    longLabel $sample coverage negative (-) strand
    type bigWig
    parent hitsclip_coverage
    color $color

    track ${posuniqbw/.bw}
    bigDataUrl $posuniqbw
    shortLabel $sample unique coverage POS
    longLabel $sample unique coverage positive (+) strand
    type bigWig
    parent hitsclip_coverage
    color $color

    track ${neguniqbw/.bw}
    bigDataUrl $neguniqbw
    shortLabel $sample unique coverage NEG
    longLabel $sample unique coverage negative (-) strand
    type bigWig
    parent hitsclip_coverage
    color $color

coverage_track
done
