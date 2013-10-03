#!/usr/bin/env bash
#BSUB -J annotatepeaks
#BSUB -e annotatepeaks.%J.err
#BSUB -o annotatepeaks.%J.out
#BSUB -q normal
#BSUB -R "select[mem>16] rusage[mem=16] span[hosts=1]"
#BSUB -n 1
#BSUB -P spillman

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/spillman/bin/config.sh

samples=(E21PR_plus E21PR_minus)
# sample=${samples[$(($LSB_JOBINDEX - 1))]}

annotatePeaks.pl E21PR_plus/peaks.txt hg19 -size 400 -d E21PR_plus E21PR_minus -m E21PR_plus_motif/all.motif > E21PR_plus/peaks.diff.annotated.txt
annotatePeaks.pl E21PR_minus/peaks.txt hg19 -size 400 -d E21PR_plus E21PR_minus -m E21PR_minus_motif/all.motif > E21PR_minus/peaks.diff.annotated.txt
annotatePeaks.pl E21PR_plus/peaks.diff.txt hg19 -size 400 -d E21PR_plus E21PR_minus -m E21PR_plus_motif/all.motif > E21PR_plus/peaks.diff.annotated.txt
annotatePeaks.pl E21PR_minus/peaks.diff.txt hg19 -size 400 -d E21PR_plus E21PR_minus -m E21PR_minus_motif/all.motif > E21PR_minus/peaks.diff.annotated.txt
