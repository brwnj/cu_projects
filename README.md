==============================================================================
Duval
==============================================================================

* UCSC Tracks
* Write email explaining steps and point to fastq, deseq plots, etc

==============================================================================
Hits-Clip
==============================================================================

* stats for snoRNA
* ayb with shorter and shorter (less priority)

==============================================================================
Marrack
==============================================================================

The ISO is an off-target IP control and Tbet is the target IP that you want to examine
input sample as a random background
then use the ISO sample to subtract any non-specific pull-downs that are seen in both ISO and Tbet

* Peaks
    * control background = RS_input_CCGTCC_L005_R1_001
* filter peaks by subtracting iso from tbet
* deliver tbet only peaks

track type=bigWig name='ISO Pos' description='ISO Coverage on positive strand' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/marrack/20121109/RS_iso_ATGTCA_L005_R1_001.pos.bw maxHeightPixels=15:50:35 color=37,52,148 visibility=full
track type=bigWig name='ISO Neg' description='ISO Coverage on negative strand' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/marrack/20121109/RS_iso_ATGTCA_L005_R1_001.neg.bw maxHeightPixels=15:50:35 color=37,52,148 visibility=full
track type=bigWig name='TBET Pos' description='TBET Coverage on positive strand' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/marrack/20121109/RS_tbet_CTTGTA_L005_R1_001.pos.bw maxHeightPixels=15:50:35 color=228,26,28 visibility=full
track type=bigWig name='TBET Neg' description='TBET Coverage on negative strand' bigDataUrl=http://amc-sandbox.ucdenver.edu/~brownj/marrack/20121109/RS_tbet_CTTGTA_L005_R1_001.neg.bw maxHeightPixels=15:50:35 color=228,26,28 visibility=full
http://amc-sandbox.ucdenver.edu/~brownj/marrack/20121109/RS_tbet_CTTGTA_L005_R1_001_unique_peaks.bed.gz


==============================================================================
Poly(A)
==============================================================================

    samtools view bam | python dexseq_count.py flattened.gff - out.counts
    awk 'BEGIN{OFS=FS="\t"}{if($2!=0){n=int(rand()*5);print $1,$2+n}else{print $0}}' out.counts > outx.counts

* neg and pos split everything; annotation, bam, etc
* tracks have to be shown only using 5' read ends
* re-run dexseq