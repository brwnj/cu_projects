#!/usr/bin/env bash
#BSUB -J intersect
#BSUB -e intersect.%J.err
#BSUB -o intersect.%J.out
#BSUB -q normal
#BSUB -R "select[mem>4] rusage[mem=4] span[hosts=1]"
#BSUB -n 1
#BSUB -P zhao

<<DOC
process intersection of sups minus wt
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/zhao/bin/config.sh

bedtools subtract -A -a sup6 -b parent > sup6.vcf
bedtools subtract -A -a sup8 -b parent > sup8.vcf

java -jar ~/opt/snpeff/snpEff.jar eff -chr chr -noStats -v -c ~/opt/snpeff/snpEff.config -o txt sacCer3 $RESULTS/intersection_minus_wt.vcf > $RESULTS/intersection_minus_wt.txt.stupidheader
tail -n +3 $RESULTS/intersection_minus_wt.txt.stupidheader > $RESULTS/intersection_minus_wt.txt
rm $RESULTS/intersection_minus_wt.txt.stupidheader
python $BIN/annotate.py $RESULTS/intersection_minus_wt.txt $ANNOTATION > $RESULTS/intersection_minus_wt.annotated.txt
