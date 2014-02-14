#!/usr/bin/env bash
#BSUB -J trackhub
#BSUB -e trackhub.%J.err
#BSUB -o trackhub.%J.out
#BSUB -q normal
#BSUB -n 1
#BSUB -P duval_exome

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/duval_exome/bin/config.sh

genome=canFam3

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
    echo "hub duval_exome" > $hub
    echo "shortLabel Duval" >> $hub
    echo "longLabel Duval" >> $hub
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
track duval_coverage
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
    color=${HUBCOLORS[$sample]}
    cp $RESULTS/$sample/*.bw $HUB/$genome
    posbw=${sample}_pos.bw
    negbw=${sample}_neg.bw
    cat <<coverage_track >>$trackdb
    track ${posbw/.bw}
    bigDataUrl $posbw
    shortLabel $sample coverage POS
    longLabel $sample coverage positive (+) strand
    type bigWig
    parent duval_coverage
    color $color

    track ${negbw/.bw}
    bigDataUrl $negbw
    shortLabel $sample coverage NEG
    longLabel $sample coverage negative (-) strand
    type bigWig
    parent duval_coverage
    color $color

coverage_track
done

# output for alignments
cat <<alignments_track >>$trackdb
track duval_alignments
compositeTrack on
shortLabel Alignments
longLabel Alignments
type bam

alignments_track

for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    cp $RESULTS/$sample/$sample.bam* $HUB/$genome
    bam=$sample.bam

    cat <<alignments_track >>$trackdb
    track ${bam/.bw}
    bigDataUrl $bam
    shortLabel $sample alignments
    longLabel $sample alignments
    type bam
    showNames off
    bamColorMode strand
    parent duval_alignments

alignments_track
done

# output for variants
cat <<variants_track >>$trackdb
track duval_variants
compositeTrack on
shortLabel Variants
longLabel Variants
type vcfTabix

variants_track

for (( i = 0; i < ${#SAMPLES[@]}; i++ )); do
    sample=${SAMPLES[$i]}
    cp $RESULTS/$sample/$sample.vcf.gz* $HUB/$genome
    vcf=$sample.vcf.gz

    cat <<variants_track >>$trackdb
    track ${vcf/.vcf.gz}
    bigDataUrl $vcf
    shortLabel $sample variants
    longLabel $sample variants
    type vcfTabix
    parent duval_variants

variants_track
done

# exome capture regions
cat <<misc_track >>$trackdb
track duval_misc
compositeTrack on
shortLabel Misc
longLabel Misc
type bigBed 3 .

misc_track

cp $HOME/ref/canFam3/canFam3.SureDesign.bb $HUB/$genome

cat <<misc_track >>$trackdb
track canFam3.SureDesign
bigDataUrl canFam3.SureDesign.bb
shortLabel canFam3 SureDesign
longLabel canFam3 SureDesign
type bigBed 6 .
color 31,120,180
thickDrawItem on
parent duval_variants

misc_track
