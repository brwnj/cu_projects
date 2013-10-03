#! /usr/bin/env bash
#BSUB -J peaks[1-6]
#BSUB -e peaks.%J.%I.err
#BSUB -o peaks.%J.%I.out
#BSUB -q short
#BSUB -R "select[mem>8] rusage[mem=8] span[hosts=1]"
#BSUB -n 1

<<DOC
call peaks using macs2.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/polya/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}
result=$RESULT/$sample
bam=$result/$sample.bam
negbam=$result/$sample.neg.bam
posbam=$result/$sample.pos.bam

negout=$result/$sample.neg
posout=$result/$sample.pos

# the output from macs2
negpeak=${negout}_peaks.bed
pospeak=${posout}_peaks.bed
negnarrowpeak=${negout}_peaks.narrowPeak
posnarrowpeak=${posout}_peaks.narrowPeak
negsummit=${negout}_summits.bed
possummit=${posout}_summits.bed
negxls=${negout}_peaks.xls
posxls=${posout}_peaks.xls

# some peaks were extending outside of genomic coords
negclipped_peak=${negout}_peaks.bed.clipped
posclipped_peak=${posout}_peaks.bed.clipped

# combined peaks with appropriate strand column
peak=$result/${sample}_peaks.bed.gz

# classifying the called peaks
posbg=$result/$sample.pos.bedgraph.gz
negbg=$result/$sample.neg.bedgraph.gz
classified=$result/${sample}_peaks.classified.bed.gz

if [[ ! -f $bam ]]; then
    echo "bam not found for $sample"
    exit
fi

# strand out the bam to facilitate strand specific peak calling
if [[ ! -f $negbam ]]; then
    samtools view -hb -f 0x10 $bam > $negbam
fi
if [[ ! -f $posbam ]]; then
    samtools view -hb -F 0x10 $bam > $posbam
fi

# call peaks on negative strand reads
if [[ ! -f $negpeak.gz ]]; then
    macs2 callpeak -t $negbam -n $negout --keep-dup auto \
        --nomodel -s 25 --extsize 5 --call-summits
    cut -f1-4 $negnarrowpeak > $negpeak
    bedClip $negpeak $CHROM_SIZES $negclipped_peak
    mv $negclipped_peak $negpeak
    gzip -f $negpeak $negnarrowpeak
    rm -f $negxls $negsummit
fi

# call peaks on positive strand reads
if [[ ! -f $pospeak.gz ]]; then
    macs2 callpeak -t $posbam -n $posout --keep-dup auto \
        --nomodel -s 25 --extsize 5 --call-summits
    cut -f1-4 $posnarrowpeak > $pospeak
    bedClip $pospeak $CHROM_SIZES $posclipped_peak
    mv $posclipped_peak $pospeak
    gzip -f $pospeak $posnarrowpeak
    rm -f $posxls $possummit
fi

# combine pos and neg, fix peak name, and fix peak strand
if [[ ! -f $peak ]]; then
    zcat $negpeak.gz $pospeak.gz \
        | awk 'BEGIN{OFS=FS="\t"}{
            split($4, basename, "/");
            $4 = basename[length(basename)];
            if($4~"neg"){$6="+"}
            if($4~"pos"){$6="-"}print}' \
        | bedtools sort -i - \
        | gzip -c > $peak
fi

# classify the peaks
if [[ ! -f $classified ]]; then
    # XXX updated the canonical region coordinates to capture missing
    # canonical sites e.g. MYC PAS needs to be class 1
    # sort to be safe
    python $BIN/classify_peaks.py --canonical-region -5 -40 $peak $posbg $negbg \
        $FASTA $CHROM_SIZES \
        | bedtools sort -i - \
        | gzip -c > $classified
fi
