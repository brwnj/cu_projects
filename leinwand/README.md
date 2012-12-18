==============================================================================
Leinwand; Mouse; miRNA
==============================================================================

    for f in *001; do echo "peaktools-identify-peaks -t ${f}_pos -s + -w 50 --trim-peaks -v genomedata > $f/${f}_pos_peaks.bed" | bsez mirna_peaks.pos; done
    for f in *001; do echo "peaktools-identify-peaks -t ${f}_neg -s - -w 50 --trim-peaks -v genomedata > $f/${f}_neg_peaks.bed" | bsez mirna_peaks.neg; done
    for f in *001; do echo "peaktools-identify-peaks -t ${f}_pos -s + -w 50 --trim-peaks --shuffle-data -v genomedata > $f/${f}_pos_shuffle_peaks.bed" | bsez mirna_peaks.pos; done
    for f in *001; do echo "peaktools-identify-peaks -t ${f}_neg -s - -w 50 --trim-peaks --shuffle-data -v genomedata > $f/${f}_neg_shuffle_peaks.bed" | bsez mirna_peaks.neg; done

    for f in E*; do echo "peaktools-identify-peaks -t $f.neg -w 50 -s - --trim-peaks -v genomedata > $f/$f_gd_neg_peaks.bed" | bsez mirna_peaks.neg; done

    /vol1/home/brownj/ref/mm9/fastas
    /vol1/home/brownj/ref/mm9/mm9.sizes
    
* peaks
* shuffled-peaks
* combine peaks - for f in *001; do cat $f/${f}_gd_shuffle_???_peaks.bed > $f/${f}_gd_shuffle_peaks.bed; done
* sort the combined-peaks - for f in *001; do bedSort $f/${f}_peaks.bed $f/${f}_peaks.bed; done
* q-values (again, i find almost nothing passing this test)
* consensus
* coverage/counts
* deseq
* make mirbase annotation for mouse
* annotate deseq results

fix the nulls

awk '$2!="None"' pos.null neg.null > peaks.null
cat pos neg > peaks.bed