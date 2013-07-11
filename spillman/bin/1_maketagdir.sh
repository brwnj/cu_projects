#!/usr/bin/env bash
#BSUB -J maketagdir
#BSUB -e maketagdir.%J.%I.err
#BSUB -o maketagdir.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P spillman

set -o nounset -o pipefail -o errexit -x

# just ran this in ../common

echo "makeTagDirectory E21PR_plus_input E2input_plus/E2input_plus.bam" | bsez maketagdir -P spillman
echo "makeTagDirectory E21PR_minus_input E2Input_minus/E2Input_minus.bam" | bsez maketagdir -P spillman
echo "makeTagDirectory E21PR_plus E21PR1_plus/E21PR1_plus.bam E21PR2_plus/E21PR2_plus.bam E21PR3_plus/E21PR3_plus.bam " | bsez maketagdir -P spillman
echo "makeTagDirectory E21PR_minus E21PR1_minus/E21PR1_minus.bam E21PR2_minus/E21PR2_minus.bam E21PR3_minus/E21PR3_minus.bam " | bsez maketagdir -P spillman