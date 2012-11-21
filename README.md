==============================================================================
Artinger
==============================================================================

* zebrafish!!
* redo: fix macs genome size!
* peaks for both groups but eventually it may be interesting to compare between the two
* ucsc tracks + common deliverables

==============================================================================
Hits-Clip
==============================================================================

* ayb with shorter and shorter (less priority)

==============================================================================
Poly(A)
==============================================================================

    samtools view bam | python dexseq_count.py flattened.gff - out.counts
    rm *x.counts
    for f in *.counts; do awk 'BEGIN{OFS=FS="\t"}!/^_/{foo=int(rand()*10);if($2!=0){print $1,$2+foo}else{print}}' $f > $(basename $f .counts)x.counts; done

* neg and pos split everything; annotation, bam, etc
* tracks have to be shown only using 5' read ends
* re-run dexseq

==============================================================================
Walter
==============================================================================

* annotate to genbank annotation file
* sequences appended to deseq results
* append each basemean for each time series
* split peaks then form consensus
* maybe re-run on human