#!/usr/bin/env bash
#BSUB -J ayb_main
#BSUB -e ayb.%J.err
#BSUB -o ayb.%J.out
#BSUB -q normal
#BSUB -P pillai_kabos_polya

<<DOC
call bases using AYB
DOC

set -o nounset -o errexit -o pipefail -x

source $HOME/projects/polya/bin/config.sh

# leaving this because sometimes there are multiple lanes
lane=3

# 50 bp
# bs=R6I14C36
# 100 bp
bs=R6I14C86
# Expected folder structure:
# $HOME/projects/polya/data/20130305/Data/Intensities/L003/<cycles>/*.cif
cifs=$HOME/projects/polya/data/20130305
# to ensure proper tiles
cifs_dir=$cifs/Data/Intensities/L00${lane}/C1.1
fastq=$cifs/20130226_L00${lane}.fastq.gz
barcodes=$cifs/barcodes.txt

jobids=""
# you need a for loop here because some seg fault
for tile in `ls $cifs_dir | sed -rn 's/._._([0-9]+).cif/\1/p'`; do
    RUNSCRIPT=ayb.${tile}.$lane.sh
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
    
    # submit the job
    job=$(bsub < $RUNSCRIPT)
    # get the id
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    # add the job id to a space delimited string
    jobids="$jobids $jobid"
done

# wait for completion
python -m bsub $jobids
# concatenate into one fastq
cat $cifs/*.fastq | gzip -c > $fastq
# demultiplex
fastq-multx -B $barcodes -m 2 -e $fastq -o %.fq
gzip $cifs/*.fq

# qc the reads of individual samples
for f in $cifs/*.fq.gz; do
    runfastqc $f $PROJECTID $cifs
done