==============================================================================
Artinger
==============================================================================

* check that the pipeline ran to completion
  
==============================================================================
Hits-Clip
==============================================================================

* stats for snoRNA
* ayb with shorter and shorter (less priority)

==============================================================================
Poly(A)
==============================================================================

    samtools view bam | python dexseq_count.py flattened.gff - out.counts
    awk 'BEGIN{OFS=FS="\t"}{if($2!=0){n=int(rand()*5);print $1,$2+n}else{print $0}}' out.counts > outx.counts

* neg and pos split everything; annotation, bam, etc
* tracks have to be shown only using 5' read ends
* re-run dexseq

==============================================================================
Walter
==============================================================================

* test across the time series
* pull out significant
* intersect with mirbase19
* compare lists

* peaks on tuberculosis
* consensus, counts, deseq
* tests across time series
* send list of anything significant