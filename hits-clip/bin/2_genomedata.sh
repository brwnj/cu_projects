#!/usr/bin/env bash
#BSUB -J gdarch
#BSUB -e gdarch.%J.err
#BSUB -o gdarch.%J.out
#BSUB -q normal
#BSUB -P pillai_kabos_hitsclip

<<DOC
create genomedata archive.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/hits-clip/bin/config.sh

bams=$RESULTS/*/*bam
for file in $bams; do
    if [[ ! -f $file.bai ]]; then
        samtools index $file
        bsub -J indexbam -o $file.out -e $file.err -P $PI -K "samtools index $file" &
    fi
done
wait

if [[ ! -d $GENOMEDATA ]]; then
    # no need to wait for this to finish
    bam2gd.py $SIZES $FASTAS $bams -o $GENOMEDATA -p pillai_kabos_hitsclip &
fi
