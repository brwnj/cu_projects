#!/usr/bin/env bash
<<DOC
Just run this script. Don't submit to LSF.
DOC

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

for track in `genomedata-info tracknames_continuous $GENOMEDATA`; do
    echo $track
    # peaks=$RESULTS/$sample/$sample.rmd.peaks.bed.gz
    # shuffledpeaks=$RESULTS/$sample/$sample.rmd.shuffle.peaks.bed.gz

    symbol="+"
    if [[ "${track#*pos}" == "$track" ]]; then
        symbol="-"
    fi

    runscript=$track.peaks.sh
    cat <<rs >$runscript
#!/usr/bin/env bash
#BSUB -J $track.peaks[1-25]
#BSUB -e $track.peaks.%J.%I.err
#BSUB -o $track.peaks.%J.%I.out
#BSUB -q normal
#BSUB -P nicoli

source $HOME/projects/nicoli/bin/config.sh
chromosome=\${CHROMOSOMES[\$((\$LSB_JOBINDEX - 1))]}
peaktools-identify-peaks -t $track -c \$chromosome -w 30 -v -s $symbol \$GENOMEDATA | gzip -c > $track.\$chromosome.peaks.bed.gz
rs
   bsub < $runscript

    runscript=$track.shuffle.peaks.sh
    cat <<rs >$runscript
#!/usr/bin/env bash
#BSUB -J $track.shuffled_peaks[1-25]
#BSUB -e $track.shuffled_peaks.%J.%I.err
#BSUB -o $track.shuffled_peaks.%J.%I.out
#BSUB -q normal
#BSUB -P nicoli

source $HOME/projects/nicoli/bin/config.sh
chromosome=\${CHROMOSOMES[\$((\$LSB_JOBINDEX - 1))]}
peaktools-identify-peaks -t $track -c \$chromosome -w 30 -v -s $symbol --shuffle-data \$GENOMEDATA | gzip -c > $track.\$chromosome.shuffle.peaks.bed.gz
rs
   bsub < $runscript
done
