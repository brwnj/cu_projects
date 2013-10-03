#! /usr/bin/env bash
#BSUB -J bg.bw[1-69]
#BSUB -e bg.bw.%J.%I.err
#BSUB -o bg.bw.%J.%I.out
#BSUB -q short
#BSUB -P pillai_kabos_polya

<<DOC
Convert aligned BAMs to bedgraph and bigwig format. Overwrites any existing files.
DOC

set -o nounset -o pipefail -o errexit -x

source $HOME/projects/polya/bin/config.sh
sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

results=$RESULT/$sample
umibam=$results/$sample.UMIs_not_removed.bam
bam=$results/$sample.bam

# samtools view -F 0x10 = forward strand
# samtools view -f 0x10 = reverse strand

strands=(pos neg)
strandargs=("-F 0x10" "-f 0x10")

for idx in ${!strands[@]}; do

    strand=${strands[$idx]}
    strandarg=${strandargs[$idx]}

    for bamfile in $bam $umibam; do

        bambase=$(basename $bamfile .bam)
        bedgraph=$results/$bambase.$strand.bedgraph
        bigwig=$results/$bambase.$strand.bw

        if [ "$strand" == "neg" ]; then
            # this adds the length of the sequence to the coordinate
            samtools view $strandarg $bamfile \
                | awk '{print $3,$4+length($10)}' \
                | sort \
                | uniq -c \
                | awk 'BEGIN{FS=" "}{print $2,$3-2,$3-1,$1}' \
                | bedtools sort -i - \
                > $bedgraph
        else
             samtools view $strandarg $bamfile \
                | awk '{print $3,$4}' \
                | sort \
                | uniq -c \
                | awk 'BEGIN{FS=" "}{print $2,$3-1,$3,$1}' \
                | bedtools sort -i - \
                > $bedgraph
        fi    
        bedGraphToBigWig $bedgraph $CHROM_SIZES $bigwig
        gzip $bedgraph
    done
done
