#!/usr/bin/env bash
#BSUB -J ayb_main
#BSUB -e ayb.%J.err
#BSUB -o ayb.%J.out
#BSUB -q normal
#BSUB -P lance

set -o nounset -o errexit -o pipefail -x

pi=lance
lane=1
bs=R100R2I2C2
cifs=$HOME/projects/lance/data/20130613

cifs_dir=$cifs/Data/Intensities/L00${lane}/C1.1
fastq=$cifs/L00${lane}.fastq.gz
jobids=""
for tile in `ls $cifs_dir | sed -rn 's/._._([0-9]+).cif/\1/p'`; do
    runscript=ayb.$tile.$lane.sh
    echo "#!/usr/bin/env bash" > $runscript
    echo "#BSUB -J ayb_slave.$tile.$lane" >> $runscript
    echo "#BSUB -R \"span[hosts=1] select[mem>8] rusage[mem=8]\"" >> $runscript
    echo "#BSUB -e ayb.%J.$lane.err" >> $runscript
    echo "#BSUB -o ayb.%J.$lane.out" >> $runscript
    echo "#BSUB -n 8" >> $runscript
    echo "#BSUB -q normal" >> $runscript
    echo "#BSUB -P $pi" >> $runscript
    echo "AYB -b $bs -d cif -l debug -o $cifs -p 8 -i $cifs -r L${lane}T${tile} --format fastq" >> $runscript
    job=$(bsub < $runscript)
    jobid=$(echo $job | sed -rn 's/.*<([0-9]+)>.*/\1/p')
    jobids="$jobids $jobid"
done
python -m bsub $jobids
cat $cifs/*${lane}*.fastq | gzip -c > $fastq
rm $cifs/?_${lane}_*.fastq
rm $cifs/ayb*.tab
