#! /usr/bin/env bash
#BSUB -J cuffdiff[1-2]
#BSUB -e cuffdiff.%J.%I.err
#BSUB -o cuffdiff.%J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>20] rusage[mem=20] span[hosts=1]"
#BSUB -n 4

<<DOC
assemble transcripts and perform differential expression testing.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/huang_hua/bin/huang_hua.cfg

# TODO: need to figure out a better way to handle test groups
pair=${TESTS[$LSB_JOBINDEX]}
testname=${TESTNAMES[$LSB_JOBINDEX]}
outdir=$RESULTS/common/$testname
label=$(echo $testname | sed 's:\_vs\_:,:')
options="-o $outdir -p4 -g $GTF -b $FASTA -u --min-reps-for-js-test 1 -L $label"
cuffdiff $options $pair

# gfold
sample1=$(echo $pair | cut -f1 -d" ")
sample2=$(echo $pair | cut -f2 -d" ")

outfile=$RESULTS/common/$testname/$testname.gfold.diff
samtools view $sample1 | gfold count -ann $GTF -annf GTF -tag stdin -o ${sample1/.bam/.read_cnt}
samtools view $sample2 | gfold count -ann $GTF -annf GTF -tag stdin -o ${sample2/.bam/.read_cnt}
gfold diff -s1 ${sample1/.bam/} -s2 ${sample2/.bam/} -suf .read_cnt -o $outfile