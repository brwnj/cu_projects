#! /usr/bin/env bash
#BSUB -J cufflinks[1-9]
#BSUB -e cufflinks.%J.%I.err
#BSUB -o cufflinks.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -n 4

<<DOC
assemble transcripts using cufflinks.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/ribera/bin/ribera.cfg

sample=${SAMPLES[$LSB_JOBINDEX]}
outdir=$RESULTS/common/$sample
options="-o $outdir -p4 -G $GTF -b $FASTA -u"
cufflinks $options $RESULTS/common/$sample/$sample.bam