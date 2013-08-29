#!/usr/bin/env bash
#BSUB -J ayb_main
#BSUB -e ayb_main.%J.err
#BSUB -o ayb_main.%J.out
#BSUB -q normal
#BSUB -P ken_ayb

<<DOC
Expected folder structure:
$HOME/projects/polya/data/20130305/Data/Intensities/L003/<cycles>/*.cif
DOC

set -o nounset -o errexit -o pipefail -x

lane=4
date=20130827

bs=R8I4C88I1R6I5C96

cifs=$HOME/projects/ken_ayb/data/$date
cifs_dir=$cifs/Data/Intensities/L00${lane}/C1.1
afastq=$cifs/${date}_L00${lane}_R1.fastq.gz
bfastq=$cifs/${date}_L00${lane}_R2.fastq.gz

jobids=""
for tile in `ls $cifs_dir | sed -rn 's/._._([0-9]+).cif/\1/p'`; do
    RUNSCRIPT=ayb.$tile.$lane.sh
    echo "#! /usr/bin/env bash" > $RUNSCRIPT
    echo "#BSUB -J ayb_slave.${tile}.$lane" >> $RUNSCRIPT
    echo "#BSUB -R \"span[hosts=1] select[mem>8] rusage[mem=8]\"" >> $RUNSCRIPT
    echo "#BSUB -e ayb.%J.$lane.err" >> $RUNSCRIPT
    echo "#BSUB -o ayb.%J.$lane.out" >> $RUNSCRIPT
    echo "#BSUB -n 8" >> $RUNSCRIPT
    echo "#BSUB -q normal" >> $RUNSCRIPT
    echo "#BSUB -P pillai_kabos_polya" >> $RUNSCRIPT
    echo "

AYB -b $bs -d cif -l debug -o $cifs -p 8 -i $cifs -r L${lane}T${tile} --format fastq
" >> $RUNSCRIPT

    job=$(bsub < $RUNSCRIPT)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
done

python -m bsub $jobids
cat $cifs/?_${lane}_*a.fastq | gzip -c > $afastq
cat $cifs/?_${lane}_*b.fastq | gzip -c > $bfastq
rm $cifs/?_${lane}_*.fastq
rm $cifs/ayb*.tab
