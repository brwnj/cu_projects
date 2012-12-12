==============================================================================
Peaktools
==============================================================================

prints None sometimes:
    
    chrY	None	None	MDX_22_AGTTCC_L003_R1_001_pos	None	+

==============================================================================
Artinger; Zebrafish; ChIP-Seq
==============================================================================

==============================================================================
Hits-Clip; hg18
==============================================================================

* of peaks in 5'UTR, 3'UTR, and CDS, identify peaks unique to control and e-treated for MCF7 and BT474

* ayb on new (20121210) polya stuff
* ayb with shorter and shorter (less priority)

==============================================================================
Leinwand; Mouse
==============================================================================

* per sample means appended to the DESeq output
* send more of the deseq output and qc including pca and clustering

==============================================================================
Jacobsen; Human; ChIP-Seq
==============================================================================

YA+EtOH & YiA+ponA+EtOH
YA+P4 & YiA+ponA+P4
YA+ZK 1h & YiA+PonA+ZK 1h
YiA EtOH & Y EtOH etc

biological replicates
2,7
3,10
4,9
2,1

==============================================================================
Poly(A)
==============================================================================

    samtools view bam | python dexseq_count.py flattened.gff - out.counts
    rm *x.counts
    for f in *.counts; do awk 'BEGIN{OFS=FS="\t"}{foo=int(rand()*10);if($2!=0){print $1,$2+foo}else{print}}' $f > $(basename $f .counts)x.counts; done


==============================================================================
Walter
==============================================================================

* rerun analysis on human samples
* find where some of the sequences are that are DE in H37Rv -- they be DE in hg19
* annotate to genome feature, primary interest is 3' UTR