#! /usr/bin/env bash
<<DOC
convert from bam to bw, gzipping all intermediary steps.
DOC

set -o nounset -o pipefail -o errexit

function usage () {
    echo "Check your args.
         "
    echo "Usage: sh $0 <BAM> <CHROMOSOME_SIZES> <GENOME>"
    echo "Example: sh $0 MP1.bam ~/reference/hg19/hg19.sizes hg19"
    exit 1
}

if [ $# -lt 3 ]; then
    usage
fi

bam=$1
chromsizes=$2
genome=$3
sample=$(basename $bam .bam)
runscript="bamtobw.sh"

cat <<runscript >${runscript}
#BSUB -J bamtobw
#BSUB -e bamtobw.%J.err
#BSUB -o bamtobw.%J.out

bedtools bamtobed -i ${sample}.bam > ${sample}.bed
awk '\$6=="+"' ${sample}.bed > ${sample}.pos.bed
awk '\$6=="-"' ${sample}.bed > ${sample}.neg.bed
bedSort ${sample}.pos.bed ${sample}.pos.bed
bedSort ${sample}.neg.bed ${sample}.neg.bed
bedItemOverlapCount -chromSize=${chromsizes} ${genome} stdin < ${sample}.pos.bed > ${sample}.pos.bedgraph
bedItemOverlapCount -chromSize=${chromsizes} ${genome} stdin < ${sample}.neg.bed > ${sample}.neg.bedgraph
bedGraphToBigWig ${sample}.pos.bedgraph ${chromsizes} ${sample}.pos.bw
bedGraphToBigWig ${sample}.neg.bedgraph ${chromsizes} ${sample}.neg.bw
gzip -f *.bed *.bedgraph
runscript

bsub < ${runscript}