#!/usr/bin/env bash
#BSUB -J gfold[1-10]
#BSUB -e %J.%I.err
#BSUB -o %J.%I.out
#BSUB -q normal
#BSUB -R "select[mem>24] rusage[mem=24] span[hosts=1]"
#BSUB -n 1

<<DOC
differential expression using gfold

[ cdiff only [ intersection ] gfold only ]
11vs12  [ 1584 [ 282 ] 98 ]
13vs14  [ 1745 [ 512 ] 210 ]
1vs2    [ 1358 [ 364 ] 33 ]
then i added significant == 'yes' criteria to filter cuffdiff output
11vs12  [ 21 [ 123 ] 257 ]
13vs14  [ 15 [ 173 ] 553 ]
1vs2    [ 1 [ 105 ] 292 ]
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/cambier/bin/cambier.cfg

pair=${TESTS[$LSB_JOBINDEX]}
sample1=$(echo $pair | cut -f1 -d" ")
sample2=$(echo $pair | cut -f2 -d" ")
testname=${TESTNAMES[$LSB_JOBINDEX]}
outdir=$RESULTS/common/$testname.diff

samtools view $sample1 | gfold count -ann $GTF -annf GTF -tag stdin -o ${sample1/.bam/.read_cnt}
samtools view $sample2 | gfold count -ann $GTF -annf GTF -tag stdin -o ${sample2/.bam/.read_cnt}
gfold diff -s1 ${sample1/.bam/} -s2 ${sample2/.bam/} -suf .read_cnt -o $outdir