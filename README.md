==============================================================================
Artinger; Zebrafish; ChIP-Seq
==============================================================================

* zebrafish!!
* redo: fix macs genome size!
* peaks for both groups but eventually it may be interesting to compare between the two
* ucsc tracks + common deliverables

==============================================================================
Hits-Clip
==============================================================================

* correlation of samples among cases
* ayb with shorter and shorter (less priority)

==============================================================================
Leinwand; Mouse; miRNA
==============================================================================

    for f in E*; do echo "peaktools-identify-peaks -t $f.neg -w 50 -s - --trim-peaks -v genomedata > $f/$f_gd_neg_peaks.bed" | bsez mirna_peaks.neg; done

* genomedata archive into pipeline. need bedgraphs for each bam first.
* need bw as well for ucsc

/vol1/home/brownj/ref/mm9/fastas
/vol1/home/brownj/ref/mm9/mm9.sizes

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

* genomedata archive!

    chr1	549717	549726	E1T1_Inf.pos	19.9870524164	+
    chr1	551813	551823	E1T1_Inf.pos	21.7511656683	+
    chr1	552338	552340	E1T1_Inf.pos	9.12588568711	+
    chr1	553000	553002	E1T1_Inf.pos	9.12588568711	+
    chr1	553106	553115	E1T1_Inf.pos	19.9870524164	+
    chr1	553478	553489	E1T1_Inf.pos	23.5609476419	+
    chr1	553850	553858	E1T1_Inf.pos	18.2707417598	+
    chr1	554500	554509	E1T1_Inf.pos	19.9870524164	+
    chr1:545523-545821
    chr1	554629	554640	E1T1_Inf.pos	23.5609476419	+
    chr1	554972	554990	E1T1_Inf.pos	37.3616240056	+
    chr1	555650	555664	E1T1_Inf.pos	29.2455161535	+
    chr1	556651	556664	E1T1_Inf.pos	27.3098517404	+
    chr1	559206	559217	E1T1_Inf.pos	23.5609476419	+
    chr1	559605	559615	E1T1_Inf.pos	21.7511656683	+

take the peaks from deseq
loop over coords, find summit location
add 30 bp up and downstream from summit
write out new best guess in addition to deseq input and sequence for the guess

* annotate to genbank annotation file
* sequences appended to deseq results
* append each basemean for each time series
* split peaks then form consensus
* maybe re-run on human