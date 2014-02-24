#!/usr/bin/env bash
#BSUB -J postprocessing
#BSUB -e postprocessing.%J.err
#BSUB -o postprocessing.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

<<DOC

DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

# sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

# combine replicate tags into single
<<tags
echo "makeTagDirectory 3B5_hela 3B5_hela_*/*.bam -genome hg19 -normGC default" | bsez -J maketagdir -P williams_chipseq
echo "makeTagDirectory 3B5_mitotic_hela 3B5_mitotic_hela_?/*.bam -genome hg19 -normGC default" | bsez -J maketagdir -P williams_chipseq
echo "makeTagDirectory polII_hela polII_hela_?/*.bam -genome hg19 -normGC default" | bsez -J maketagdir -P williams_chipseq
echo "makeTagDirectory polII_mitotic_hela polII_mitotic_hela_?/*.bam -genome hg19 -normGC default" | bsez -J maketagdir -P williams_chipseq
echo "makeTagDirectory SC184_hela SC184_hela_?/*.bam -genome hg19 -normGC default" | bsez -J maketagdir -P williams_chipseq
echo "makeTagDirectory SC184_mitotic_hela SC184_mitotic_hela_?/*.bam -genome hg19 -normGC default" | bsez -J maketagdir -P williams_chipseq
tags

<<peaks
echo "findPeaks 3B5_hela -style factor -o auto -i hela_input" | bsez -J findpeaks -P williams_chipseq
echo "findPeaks 3B5_mitotic_hela -style factor -o auto -i mitotic_hela_input" | bsez -J findpeaks -P williams_chipseq
echo "findPeaks polII_hela -style factor -o auto -i hela_input" | bsez -J findpeaks -P williams_chipseq
echo "findPeaks polII_mitotic_hela -style factor -o auto -i mitotic_hela_input" | bsez -J findpeaks -P williams_chipseq
echo "findPeaks SC184_hela -style factor -o auto -i hela_input" | bsez -J findpeaks -P williams_chipseq
echo "findPeaks SC184_mitotic_hela -style factor -o auto -i mitotic_hela_input" | bsez -J findpeaks -P williams_chipseq
peaks

# mergePeaks -prefix merge 3B5_hela/peaks.txt SC184_hela/peaks.txt
# mergePeaks -prefix merge 3B5_mitotic_hela/peaks.txt SC184_mitotic_hela/peaks.txt
# mergePeaks -prefix merge 3B5_hela/peaks.txt 3B/peaks.txt
# mergePeaks -prefix merge 3B5_hela/peaks.txt 3B5_mitotic_hela/peaks.txt
# mergePeaks -prefix merge SC184_hela/peaks.txt SC184_mitotic_hela/peaks.txt

# differentially bound
echo "getDifferentialPeaks merge_3B5_hela_peaks.txt_SC184_hela_peaks.txt 3B5_hela SC184_hela -F 2 -P 0.0001 > 3B5_hela__SC184_hela_diff.txt" | bsez -J diffboundpeaks -P williams_chipseq
echo "getDifferentialPeaks merge_3B5_hela_peaks.txt_SC184_hela_peaks.txt 3B5_hela SC184_hela -F 2 -P 0.0001 -same > 3B5_hela__SC184_hela_same.txt" | bsez -J diffboundpeaks -P williams_chipseq
echo "getDifferentialPeaks merge_3B5_mitotic_hela_peaks.txt_SC184_mitotic_hela_peaks.txt 3B5_mitotic_hela SC184_mitotic_hela -F 2 -P 0.0001 > 3B5_mitotic_hela__SC184_mitotic_hela_diff.txt" | bsez -J diffboundpeaks -P williams_chipseq
echo "getDifferentialPeaks merge_3B5_mitotic_hela_peaks.txt_SC184_mitotic_hela_peaks.txt 3B5_mitotic_hela SC184_mitotic_hela -F 2 -P 0.0001 -same > 3B5_mitotic_hela__SC184_mitotic_hela_same.txt" | bsez -J diffboundpeaks -P williams_chipseq
echo "getDifferentialPeaks merge_3B5_hela_peaks.txt_3B5_mitotic_hela_peaks.txt 3B5_hela 3B5_mitotic_hela -F 2 -P 0.0001 > 3B5_hela__3B5_mitotic_hela_diff.txt" | bsez -J diffboundpeaks -P williams_chipseq
echo "getDifferentialPeaks merge_3B5_hela_peaks.txt_3B5_mitotic_hela_peaks.txt 3B5_hela 3B5_mitotic_hela -F 2 -P 0.0001 -same > 3B5_hela__3B5_mitotic_hela_same.txt" | bsez -J diffboundpeaks -P williams_chipseq
echo "getDifferentialPeaks merge_SC184_hela_peaks.txt_SC184_mitotic_hela_peaks.txt SC184_hela SC184_mitotic_hela -F 2 -P 0.0001 > SC184_hela__SC184_mitotic_hela_diff.txt" | bsez -J diffboundpeaks -P williams_chipseq
echo "getDifferentialPeaks merge_SC184_hela_peaks.txt_SC184_mitotic_hela_peaks.txt SC184_hela SC184_mitotic_hela -F 2 -P 0.0001 -same > SC184_hela__SC184_mitotic_hela_same.txt" | bsez -J diffboundpeaks -P williams_chipseq

# annotate diff bound peaks
echo "annotatePeaks.pl 3B5_hela__SC184_hela_diff.txt hg19 -d 3B5_hela SC184_hela > 3B5_hela__SC184_hela_diff_annotated.txt" | bsez -P williams_chipseq -J annotatepeaks
echo "annotatePeaks.pl 3B5_hela__SC184_hela_same.txt hg19 -d 3B5_hela SC184_hela > 3B5_hela__SC184_hela_same_annotated.txt" | bsez -P williams_chipseq -J annotatepeaks
echo "annotatePeaks.pl 3B5_mitotic_hela__SC184_mitotic_hela_diff.txt hg19 -d 3B5_mitotic_hela SC184_mitotic_hela > 3B5_mitotic_hela__SC184_mitotic_hela_diff_annotated.txt" | bsez -P williams_chipseq -J annotatepeaks
echo "annotatePeaks.pl 3B5_mitotic_hela__SC184_mitotic_hela_same.txt hg19 -d 3B5_mitotic_hela SC184_mitotic_hela > 3B5_mitotic_hela__SC184_mitotic_hela_same_annotated.txt" | bsez -P williams_chipseq -J annotatepeaks
echo "annotatePeaks.pl 3B5_hela__3B5_mitotic_hela_diff.txt hg19 -d 3B5_hela 3B5_mitotic_hela > 3B5_hela__3B5_mitotic_hela_diff_annotated.txt" | bsez -P williams_chipseq -J annotatepeaks
echo "annotatePeaks.pl 3B5_hela__3B5_mitotic_hela_same.txt hg19 -d 3B5_hela 3B5_mitotic_hela > 3B5_hela__3B5_mitotic_hela_same_annotated.txt" | bsez -P williams_chipseq -J annotatepeaks
echo "annotatePeaks.pl SC184_hela__SC184_mitotic_hela_diff.txt hg19 -d SC184_hela SC184_mitotic_hela > SC184_hela__SC184_mitotic_hela_diff_annotated.txt" | bsez -P williams_chipseq -J annotatepeaks
echo "annotatePeaks.pl SC184_hela__SC184_mitotic_hela_same.txt hg19 -d SC184_hela SC184_mitotic_hela > SC184_hela__SC184_mitotic_hela_same_annotated.txt" | bsez -P williams_chipseq -J annotatepeaks
