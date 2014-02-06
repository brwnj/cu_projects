#!/usr/bin/env bash
#BSUB -J annotatepeaks[1-20]
#BSUB -e annotatepeaks.%J.%I.err
#BSUB -o annotatepeaks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P williams_chipseq

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/williams_chipseq/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

if [[ $sample == *input* ]]; then
    exit
fi

peaks=$RESULTS/$sample/peaks.txt
motifsresultfile=$RESULTS/$sample/knownResults.txt
motifs=$RESULTS/$sample/known_motifs.txt
annotatedpeaks=$RESULTS/$sample/annotated_peaks.txt
replicate=${sample%??}

if [[ ! -f $motifsresultfile ]]; then
    exit
fi

if [[ ! -f $motifs ]]; then
    awk -t 'NR>1{print $2}' $motifsresultfile > $motifs
fi

# for each individual sample, use its peaks file, but always compare it to the
# other 2 reps
rep1=$RESULTS/${replicate}_1
rep2=$RESULTS/${replicate}_2
rep3=$RESULTS/${replicate}_3

annotatePeaks.pl $peaks hg19 -d $rep1 $rep2 $rep3 -m $motifs > $annotatedpeaks
