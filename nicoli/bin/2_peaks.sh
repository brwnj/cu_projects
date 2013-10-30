#!/usr/bin/env bash
#BSUB -J callpeaks[1-14]
#BSUB -e callpeaks.%J.%I.err
#BSUB -o callpeaks.%J.%I.out
#BSUB -q normal
#BSUB -P nicoli

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

peaks=$RESULTS/$sample/$sample.rmd.peaks.bed.gz
shuffledpeaks=$RESULTS/$sample/$sample.rmd.shuffle.peaks.bed.gz

if [ -f $peaks ] && [ -f $shuffledpeaks ]; then
    echo "processing complete for $sample"
    exit 1
fi

for strand in pos neg; do
    symbol="+"
    if [[ "$strand" == "neg" ]]; then
        symbol="-"
    fi
    # regular peaks
    runscript="$sample.$strand.peaks.sh"
    echo "#!/usr/bin/env bash" > $runscript
    echo "#BSUB -J $sample.$strand.peaks[1-25]" >> $runscript
    echo "#BSUB -e $sample.$strand.%J.%I.err" >> $runscript
    echo "#BSUB -o $sample.$strand.%J.%I.out" >> $runscript
    echo "#BSUB -q normal" >> $runscript
    echo "#BSUB -P nicoli" >> $runscript
    echo "
chromosomes=($CHROMOSOMES)
chromosome=\${chromosomes[\$((\$LSB_JOBINDEX - 1))]}
peaktools-identify-peaks -t ${sample}.rmd_${strand} -c \$chromosome -w 30 -v -s $symbol $GENOMEDATA | gzip -c > $sample.\$chromosome.$strand.rmd.peaks.bed.gz
" >> $runscript
    bsub < $runscript

    #shuffled peaks
    runscript="$sample.$strand.shuffle.peaks.sh"
    echo "#! /usr/bin/env bash" > $runscript
    echo "#BSUB -J $sample.$strand.shuffle[1-25]" >> $runscript
    echo "#BSUB -e %J.%I.err" >> $runscript
    echo "#BSUB -o %J.%I.out" >> $runscript
    echo "#BSUB -q normal" >> $runscript
    echo "#BSUB -P nicoli" >> $runscript
    echo "
chromosomes=($CHROMOSOMES)
chromosome=\${chromosomes[\$((\$LSB_JOBINDEX - 1))]}
peaktools-identify-peaks -t ${sample}.rmd_${strand} -c \$chromosome -w 30 -v -s $symbol --shuffle-data $GENOMEDATA | gzip -c > $sample.\$chromosome.$strand.rmd.shuffle.peaks.bed.gz
" >> $runscript
    bsub < $runscript
done
