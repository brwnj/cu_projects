#!/usr/bin/env bash
#BSUB -J callpeaks[1-10]
#BSUB -e callpeaks.%J.%I.err
#BSUB -o callpeaks.%J.%I.out
#BSUB -q normal
#BSUB -P nicoli

set -o nounset -o pipefail -o errexit -x
source $HOME/projects/nicoli/bin/config.sh

sample=${SAMPLES[$(($LSB_JOBINDEX - 1))]}

for strand in pos neg; do
    for symbol in + -; do
        # regular peaks
        runscript="$sample.$strand.peaks.sh"
        echo "#!/usr/bin/env bash" > $runscript
        echo "#BSUB -J $sample.$strand.peaks[1-25]" >> $runscript
        echo "#BSUB -e $sample.$strand.%J.%I.err" >> $runscript
        echo "#BSUB -o $sample.$strand.%J.%I.out" >> $runscript
        echo "#BSUB -q normal" >> $runscript
        echo "#BSUB -P nicoli" >> $runscript
        echo "
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
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
chromosomes=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM)
chromosome=\${chromosomes[\$((\$LSB_JOBINDEX - 1))]}
peaktools-identify-peaks -t ${sample}.rmd_${strand} -c \$chromosome -w 30 -v -s $symbol --shuffle-data $GENOMEDATA | gzip -c > $sample.\$chromosome.$strand.rmd.shuffle.peaks.bed.gz
" >> $runscript
        bsub < $runscript
    done
done
